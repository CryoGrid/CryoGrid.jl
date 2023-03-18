const RESERVED_LAYER_NAMES = (:top, :bottom, :strat, :init, :event)

# NamedLayer type (alias of Named for Layer types)
const NamedLayer{name,TLayer} = Named{name,TLayer} where {name,TLayer<:Layer}
@inline layertype(layer::NamedLayer) = layertype(typeof(layer))
@inline layertype(::Type{<:NamedLayer{name,TLayer}}) where {name,TLayer} = TLayer
@inline layername(layer::NamedLayer) = layername(typeof(layer))
@inline layername(::Type{<:NamedLayer{name}}) where {name} = name

"""
    Stratigraphy{N,TLayers,TBoundaries}

Defines a 1-dimensional stratigraphy by connecting a top and bottom layer to one or more subsurface layers.
"""
struct Stratigraphy{N,TLayers,TBoundaries}
    boundaries::TBoundaries
    layers::TLayers
    Stratigraphy(boundaries::NTuple{N,Any}, layers::NTuple{N,NamedLayer}) where {N} = new{N,typeof(layers),typeof(boundaries)}(boundaries, layers)
    Stratigraphy(
        top::Pair{<:DistQuantity,<:Top},
        sub::Pair{<:DistQuantity,<:Pair{Symbol,<:SubSurface}},
        bot::Pair{<:DistQuantity,<:Bottom}
    ) = Stratigraphy(top,(sub,),bot)
    Stratigraphy(
        # use @nospecialize to (hopefully) reduce compilation overhead
        @nospecialize(top::Pair{<:DistQuantity,<:Top}),
        @nospecialize(sub::AbstractVector{<:Pair{<:DistQuantity,<:Pair{Symbol,<:SubSurface}}}),
        @nospecialize(bot::Pair{<:DistQuantity,<:Bottom})
    ) = Stratigraphy(top, Tuple(sub), bot)
    function Stratigraphy(
        # use @nospecialize to (hopefully) reduce compilation overhead
        @nospecialize(top::Pair{<:DistQuantity,<:Top}),
        @nospecialize(sub::Tuple{Vararg{Pair{<:DistQuantity,<:Pair{Symbol,<:SubSurface}}}}),
        @nospecialize(bot::Pair{<:DistQuantity,<:Bottom})
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
layers(strat::Stratigraphy) = getfield(strat, :layers)
boundaries(strat::Stratigraphy) = getfield(strat, :boundaries)
boundarypairs(strat::Stratigraphy, z_bottom) = boundarypairs(boundaries(strat), z_bottom)
boundarypairs(bounds::NTuple, z_bottom) = tuplejoin(map(tuple, bounds[1:end-1], bounds[2:end]), ((bounds[end], z_bottom),))
layernames(strat::Stratigraphy) = map(layername, layers(strat))
layertypes(::Type{<:Stratigraphy{N,TLayers}}) where {N,TLayers} = map(layertype, TLayers.parameters)
Base.keys(strat::Stratigraphy) = layernames(strat)
Base.values(strat::Stratigraphy) = layers(strat)
Base.propertynames(strat::Stratigraphy) = Base.keys(strat)
Base.getproperty(strat::Stratigraphy, sym::Symbol) = strat[Val{sym}()].val
Base.getindex(strat::Stratigraphy, sym::Symbol) = strat[Val{sym}()].val
@generated Base.getindex(strat::Stratigraphy{N,TC}, ::Val{sym}) where {N,TC,sym} = :(layers(strat)[$(findfirst(T -> layername(T) == sym, TC.parameters))])
# Array and iteration overrides
Base.size(strat::Stratigraphy) = size(layers(strat))
Base.length(strat::Stratigraphy) = length(layers(strat))
Base.getindex(strat::Stratigraphy, i::Int) = layers(strat)[i]
Base.iterate(strat::Stratigraphy) = (layers(strat)[1],layers(strat)[2:end])
Base.iterate(strat::Stratigraphy, itrstate::Tuple) = (itrstate[1],itrstate[2:end])
Base.iterate(strat::Stratigraphy, itrstate::Tuple{}) = nothing
Base.show(io::IO, strat::Stratigraphy) = print(io, "Stratigraphy($(prod("$b => $n" for (n,b) in zip(layers(strat), boundaries(strat))))")
Base.NamedTuple(strat::Stratigraphy) = (; map(named_layer -> nameof(named_layer) => named_layer, layers(strat))...)
# ConstructionBase
ConstructionBase.getproperties(strat::Stratigraphy) = (;map(Pair, Base.keys(strat), Base.values(strat))...)
function ConstructionBase.setproperties(strat::Stratigraphy, patch::NamedTuple)
    layers_patched = map(layers(strat)) do layer
        Named(layername(layer), get(patch, layername(layer), layer.val))
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
    can_interact_exprs = map(i -> :(CryoGrid.isactive(strat[$i].val, getproperty(state, layername(strat[$i])))), tuple(1:N...))
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
                    f!(strat[$i].val, strat[$j].val, getproperty(state, names[$i]), getproperty(state, names[$j]))
                end
            else
                innerchecks_exprs = j > i+1 ? map(k -> :(can_interact[$k]), i+1:j-1) : :(false)
                quote
                    if can_interact[$i] && can_interact[$j] && !any(tuple($(innerchecks_exprs...)))
                        f!(strat[$i].val, strat[$j].val, getproperty(state, names[$i]), getproperty(state, names[$j]))
                    end
                end
            end
            push!(expr.args, expr_next)
        end
    end
    push!(expr.args, :(return nothing))
    return expr
end

@inline function CryoGrid.initialcondition!(strat::Stratigraphy, state::TileState, inits)
    # initialcondition! is only called once so we don't need to worry about performance;
    # we can just loop over everything naively
    for i in 1:length(strat)-1
        layerᵢ = strat[i].val
        layerᵢ₊₁ = strat[i+1].val
        stateᵢ = getproperty(state, layername(strat[i]))
        stateᵢ₊₁ = getproperty(state, layername(strat[i+1]))
        # first invoke initialcondition! with initializers
        for init in inits
            if i == 1 && haskey(stateᵢ.states, varname(init))
                CryoGrid.initialcondition!(layerᵢ, stateᵢ, init)
            end
            if haskey(stateᵢ₊₁.states, varname(init))
                CryoGrid.initialcondition!(layerᵢ₊₁, stateᵢ₊₁, init)
            end
            if haskey(stateᵢ.states, varname(init)) && haskey(stateᵢ₊₁.states, varname(init))
                CryoGrid.initialcondition!(layerᵢ, layerᵢ₊₁, stateᵢ, stateᵢ₊₁, init)
            end
        end
        # then invoke initialcondition! standalone
        if i == 1
            CryoGrid.initialcondition!(layerᵢ, stateᵢ)
        end
        CryoGrid.initialcondition!(layerᵢ₊₁, stateᵢ₊₁)
        CryoGrid.initialcondition!(layerᵢ, layerᵢ₊₁, stateᵢ, stateᵢ₊₁)
    end
end

@inline function CryoGrid.diagnosticstep!(strat::Stratigraphy, state::TileState)
    fastiterate(layers(strat)) do named_layer
        CryoGrid.diagnosticstep!(named_layer.val, getproperty(state, layername(named_layer)))
    end
end

@inline function CryoGrid.interact!(strat::Stratigraphy, state::TileState)
    # interact! requires special implementation via `stratiterate`
    # this allows for layer states to determine which adjacent layers can and cannot interact
    stratiterate(strat, state) do layer1, layer2, state1, state2
        CryoGrid.interact!(layer1, layer2, state1, state2)
    end
end

@inline function CryoGrid.prognosticstep!(strat::Stratigraphy, state::TileState)
    fastiterate(layers(strat)) do named_layer
        CryoGrid.prognosticstep!(named_layer.val, getproperty(state, layername(named_layer)))
    end
end

@inline function CryoGrid.timestep(strat::Stratigraphy, state::TileState)
    max_dts = fastmap(layers(strat)) do named_layer
        CryoGrid.timestep(named_layer.val, getproperty(state, layername(named_layer)))
    end
    return minimum(max_dts)
end

CryoGrid.hasfixedvolume(::Type{TStrat}) where {TStrat<:Stratigraphy} = return all(map(CryoGrid.hasfixedvolume, layertypes(TStrat)))

# collecting/grouping components
CryoGrid.events(strat::Stratigraphy) = map(named_layer -> _addlayerfield(CryoGrid.events(named_layer.val), nameof(named_layer)), NamedTuple(strat))
function CryoGrid.variables(strat::Stratigraphy)
    layervars = map(CryoGrid.variables, NamedTuple(strat))
    stratvars = if !CryoGrid.hasfixedvolume(typeof(strat))
        # otherwise, if any layer has dynamic volume, make layer boundaries prognostic
        map(NamedTuple(strat)) do _
            (
                Prognostic(Symbol(:z), Scalar, u"m"),
            )
        end
    else
        # return empty tuple (no additional variables) for each layer
        map(layer -> (), NamedTuple(strat))
    end
    return map(layervars, stratvars) do lv, sv
        tuplejoin(lv, sv)
    end
end
function CryoGrid.variables(@nospecialize(named_layer::NamedLayer))
    layer = named_layer.val
    declared_vars = variables(layer)
    nested_vars = Flatten.flatten(layer, Flatten.flattenable, Var)
    all_vars = tuplejoin(declared_vars, nested_vars)
    # check for (permissible) duplicates between variables, excluding parameters
    groups = Utils.groupby(var -> varname(var), all_vars)
    for (id,vargroup) in filter(g -> length(g.second) > 1, groups)
        # if any duplicate variable deifnitions do not match, raise an error
        @assert all(vargroup[i] == vargroup[i-1] for i in 2:length(vargroup)) "Found one or more conflicting definitions of $id in $vargroup"
    end
    diag_vars = filter(isdiagnostic, all_vars)
    prog_vars = filter(isprognostic, all_vars)
    alg_vars = filter(isalgebraic, all_vars)
    # check for duplicated algebraic/prognostic vars
    prog_alg_duplicated = prog_vars ∩ alg_vars
    @assert isempty(prog_alg_duplicated) "Variables $(prog_alg_duplicated) cannot be both prognostic and algebraic."
    # check for re-definition of diagnostic variables as prognostic
    prog_alg = prog_vars ∪ alg_vars
    diag_prog = filter(v -> v ∈ prog_alg, diag_vars)
    # check for conflicting definitions of differential vars
    diff_varnames = map(v -> varname(Delta(v)), prog_alg)
    @assert all((isempty(filter(v -> varname(v) == d, all_vars)) for d in diff_varnames)) "Variable names $(Tuple(diff_varnames)) are reserved for differentials."
    # prognostic takes precedence, so we remove duplicated variables from the diagnostic variable set
    diag_vars = filter(v -> v ∉ diag_prog, diag_vars)
    # filter remaining duplicates
    diag_vars = unique(diag_vars)
    prog_vars = unique(prog_vars)
    alg_vars = unique(alg_vars)
    # convert back to tuples
    diag_vars, prog_vars, alg_vars = Tuple(diag_vars), Tuple(prog_vars), Tuple(alg_vars)
    return tuplejoin(diag_vars, prog_vars, alg_vars)
end

"""
Rebuilds the `obj` adding `name` to the `layer` field to all `Param`s, if any are defined.
"""
function _addlayerfield(@nospecialize(obj), name::Symbol)
    params = ModelParameters.params(obj)
    if length(params) > 0
        # create sub-model and add layer name to all parameters
        m = Model(obj)
        m[:layer] = repeat([name], length(params))
        return parent(m)
    else
        return obj
    end
end
