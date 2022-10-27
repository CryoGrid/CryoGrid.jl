const RESERVED_LAYER_NAMES = (:top, :bottom, :strat, :init, :event)

const NamedLayer{name,TLayer} = Named{name,TLayer} where {name,TLayer<:Layer}
@inline layertype(layer::NamedLayer) = layertype(typeof(layer))
@inline layertype(::Type{<:NamedLayer{name,TLayer}}) where {name,TLayer} = TLayer
@inline layername(layer::NamedLayer) = layername(typeof(layer))
@inline layername(::Type{<:NamedLayer{name}}) where {name} = name

"""
Type bound for stratigraphy boundaries. Boundaries may be specified as fixed distance quantities
or via some arbitrary parameterization.
"""
const StratBoundaryType = Union{<:DistQuantity,<:AbstractParam,<:Parameterization}

"""
    Stratigraphy{N,TLayers,TBoundaries}

Defines a 1-dimensional stratigraphy by connecting a top and bottom layer to one or more subsurface layers.
"""
struct Stratigraphy{N,TLayers,TBoundaries}
    boundaries::TBoundaries
    layers::TLayers
    Stratigraphy(boundaries::NTuple{N,Any}, layers::NTuple{N,NamedLayer}) where {N} = new{N,typeof(layers),typeof(boundaries)}(boundaries, layers)
    Stratigraphy(
        top::Pair{<:StratBoundaryType,<:Top},
        sub::Pair{<:StratBoundaryType,<:Pair{Symbol,<:SubSurface}},
        bot::Pair{<:StratBoundaryType,<:Bottom}
    ) = Stratigraphy(top,(sub,),bot)
    function Stratigraphy(
        # use @nospecialize to (hopefully) reduce compilation overhead
        @nospecialize(top::Pair{<:StratBoundaryType,<:Top}),
        @nospecialize(sub::Tuple{Vararg{Pair{<:StratBoundaryType,<:Pair{Symbol,<:SubSurface}}}}),
        @nospecialize(bot::Pair{<:StratBoundaryType,<:Bottom})
    )
        @assert length(sub) > 0 "At least one subsurface layer must be specified"
        top = top[1] => Named(:top => top[2])
        bot = bot[1] => Named(:bottom => bot[2])
        sub = map(x -> x[1] => Named(x[2]), sub)
        names = map(nameof, map(last, sub))
        @assert length(unique(names)) == length(names) "All layer names in Stratigraphy must be unique"
        boundaries = Tuple(map(first, (top, sub..., bot)))
        @assert issorted(boundaries) "Stratigraphy boundary locations must be in strictly increasing order."
        # wrap layers
        layers = tuple(last(top), map(last, sub)..., last(bot))
        # construct type
        new{length(layers),typeof(layers),typeof(boundaries)}(boundaries, layers)
    end
end
"""
Convenience macro for defining stratigraphies with multiple subsurface layers.
"""
macro Stratigraphy(args...)
    @assert length(args) >= 3 "At least three stratigraphy nodes (top, subsurface, bottom) must be provided!"
    if length(args) == 3
        :(Stratigraphy($(esc(args[1])), $(esc(args[2])), $(esc(args[3]))))
    elseif length(args) > 3
        :(Stratigraphy($(esc(args[1])), tuple($(esc.(args[2:end-1])...)), $(esc(args[end]))))
    end
end
@inline layers(strat::Stratigraphy) = getfield(strat, :layers)
@inline boundaries(strat::Stratigraphy) = getfield(strat, :boundaries)
@inline boundarypairs(strat::Stratigraphy, z_bottom) = boundarypairs(boundaries(strat), z_bottom)
@inline boundarypairs(bounds::NTuple, z_bottom) = tuplejoin(map(tuple, bounds[1:end-1], bounds[2:end]), ((bounds[end], z_bottom),))
@inline layernames(strat::Stratigraphy) = map(layername, layers(strat))
@inline layertypes(::Type{<:Stratigraphy{N,TLayers}}) where {N,TLayers} = map(layertype, TLayers.parameters)
@inline Base.keys(strat::Stratigraphy) = layernames(strat)
@inline Base.values(strat::Stratigraphy) = layers(strat)
@inline Base.propertynames(strat::Stratigraphy) = Base.keys(strat)
@inline Base.getproperty(strat::Stratigraphy, sym::Symbol) = strat[Val{sym}()].obj
@inline Base.getindex(strat::Stratigraphy, sym::Symbol) = strat[Val{sym}()].obj
@generated Base.getindex(strat::Stratigraphy{N,TC}, ::Val{sym}) where {N,TC,sym} = :(layers(strat)[$(findfirst(T -> layername(T) == sym, TC.parameters))])
# Array and iteration overrides
Base.size(strat::Stratigraphy) = size(layers(strat))
Base.length(strat::Stratigraphy) = length(layers(strat))
Base.getindex(strat::Stratigraphy, i::Int) = layers(strat)[i]
Base.iterate(strat::Stratigraphy) = (layers(strat)[1],layers(strat)[2:end])
Base.iterate(strat::Stratigraphy, itrstate::Tuple) = (itrstate[1],itrstate[2:end])
Base.iterate(strat::Stratigraphy, itrstate::Tuple{}) = nothing
Base.show(io::IO, strat::Stratigraphy) = print(io, "Stratigraphy($(prod("$b => $n" for (n,b) in zip(layers(strat), boundaries(strat))))")
# ConstructionBase
ConstructionBase.getproperties(strat::Stratigraphy) = (;map(Pair, Base.keys(strat), Base.values(strat))...)
function ConstructionBase.setproperties(strat::Stratigraphy, patch::NamedTuple)
    layers_patched = map(layers(strat)) do layer
        get(patch, layername(layer), layer)
    end
    return Stratigraphy(boundaries(strat), layers_patched)
end
"""
    stratiterate(f!::F, strat::Stratigraphy{N,TLayers}, state) where {F,N,TLayers}

`stratiterate` invokes the user-supplied function `f!(layer1, layer2, state1, state2)` on each pair of layers
in the stratigraphy which are adjacent and "active" based on the current `state`. `state` must have properties defined
corresponding to the name of each layer such that `getproperty(state, layername(strat[i]))` would return the appropriate
state object for the i'th layer in the stratigraphy.
"""
@generated function stratiterate(f!::F, strat::Stratigraphy{N,TLayers}, state) where {F,N,TLayers}
    expr = Expr(:block)
    # build expressions for checking whether each layer is active
    can_interact_exprs = map(i -> :(CryoGrid.thickness(strat[$i].obj, getproperty(state, layername(strat[$i]))) > 0), tuple(1:N...))
    push!(
        expr.args,
        quote
            names = layernames(strat)
            can_interact = tuple($(can_interact_exprs...))
        end
    )
    # We only invoke interact! on pairs of layers for which the following are satisfied:
    # 1) Both layers have thickness > 0 (i.e. they occupy non-zero space in the stratigraphy)
    # 2) Both layers are adjacent, or more crudely, all layers in between (if any) have zero thickness;
    # In order to make this type stable, we pre-generate all possible interact! expressions in the
    # downward direction and add a check to each one.
    for i in 1:N-1
        for j in i+1:N
            expr_next = if (i == 1 && j == 2) || (i == N-1 && j == N)
                # always apply top and bottom interactions
                quote
                    f!(strat[$i].obj, strat[$j].obj, getproperty(state, names[$i]), getproperty(state, names[$j]))
                end
            else
                innerchecks_exprs = j > i+1 ? map(k -> :(can_interact[$k]), i+1:j-1) : :(false)
                quote
                    if can_interact[$i] && can_interact[$j] && !any(tuple($(innerchecks_exprs...)))
                        f!(strat[$i].obj, strat[$j].obj, getproperty(state, names[$i]), getproperty(state, names[$j]))
                    end
                end
            end
            push!(expr.args, expr_next)
        end
    end
    push!(expr.args, :(return nothing))
    return expr
end
