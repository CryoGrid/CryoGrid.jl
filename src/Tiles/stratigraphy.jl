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
boundarypairs(strat::Stratigraphy) = boundarypairs(boundaries(strat))
boundarypairs(bounds::NTuple) = tuple(map(tuple, bounds[1:end-1], bounds[2:end])..., (bounds[end], bounds[end]))
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
    is_active_exprs = map(i -> :(CryoGrid.isactive(strat[$i].val, getproperty(state, layername(strat[$i])))), tuple(1:N...))
    # header code; get layer names and evaluate `isactive` for each layer
    push!(
        expr.args,
        quote
            names = layernames(strat)
            is_active = tuple($(is_active_exprs...))
        end
    )
    # We only invoke the function on pairs of layers for which the following are satisfied:
    # 1) Both layers are currently "active" according to CryoGrid.isactive
    # 2) Both layers are adjacent or all layers imbetween are currently inactive
    # In order to make this type stable, we generate code for all combinations of layers i, j where j > i.
    # This is effectively just a form of loop unrolling which greatly improves performance and avoids allocations.
    for i in 1:N-1
        for j in i+1:N
            expr_ij =
                if j == i + 1
                    # always invoke for adjacent layers
                    quote
                        f!(strat[$i].val, strat[$j].val, getproperty(state, names[$i]), getproperty(state, names[$j]))
                    end
                else
                    # for layers that are non-adjacent, we only invoke f! if all layers between them are currently inactive.
                    inactive_imbetween_layers = map(i+1:j-1) do k
                        quote
                            !is_active[$k]
                        end
                    end
                    quote
                        # if layers i and j are both active, and all imbetween layers are inactive, invoke f!;
                        # note the splat syntax here: $(inactive_imbetween_layers...) simply expands the tuple of
                        # expressions 'inactive_imbetween_layers' as arguments to `tuple`.
                        if is_active[$i] && is_active[$j] && all(tuple($(inactive_imbetween_layers...)))
                            f!(strat[$i].val, strat[$j].val, getproperty(state, names[$i]), getproperty(state, names[$j]))
                        end
                    end
                end
            push!(expr.args, expr_ij)
        end
    end
    push!(expr.args, :(return nothing))
    return expr
end

function Numerics.makegrid(strat::Stratigraphy, strategy::DiscretizationStrategy)
    strat_grid = nothing
    for (bounds, named_layer) in zip(boundarypairs(strat)[2:end-1], layers(strat)[2:end-1])
        if bounds[2] - bounds[1] <= zero(bounds[1])
            continue
        end
        layer = named_layer.val
        layer_grid = Numerics.makegrid(layer, strategy, bounds)
        if !isnothing(strat_grid)
            # check that grid edges line up at layer boundary
            @assert strat_grid[end] == layer_grid[1] "Upper boundary of layer $(nameof(named_layer)) does not match the previous layer."
            # concatenate grids, omitting first value of layer_grid to avoid duplicating the shared edge
            strat_grid = Grid(vcat(strat_grid, layer_grid[2:end]))
        else
            strat_grid = layer_grid
        end
    end
    return strat_grid
end

CryoGrid.initializers(strat::Stratigraphy) = tuplejoin(map(initializers, map(l -> l.val, layers(strat)))...)

function CryoGrid.initialcondition!(strat::Stratigraphy, state, inits)
    # initialcondition! is only called once so we don't need to worry about performance;
    # we can just loop over everything naively
    all_inits = tuple(initializers(strat)..., inits...)
    for i in 1:length(strat)-1
        layerᵢ = strat[i].val
        layerᵢ₊₁ = strat[i+1].val
        stateᵢ = getproperty(state, layername(strat[i]))
        stateᵢ₊₁ = getproperty(state, layername(strat[i+1]))
        # first invoke initialcondition! with initializers
        for init in all_inits
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

@inline function CryoGrid.updatestate!(strat::Stratigraphy, state)
    fastiterate(layers(strat)) do named_layer
        CryoGrid.updatestate!(named_layer.val, getproperty(state, layername(named_layer)))
    end
end

@inline function CryoGrid.interact!(strat::Stratigraphy, state)
    # interact! requires special implementation via `stratiterate`
    # this allows for layer states to determine which adjacent layers can and cannot interact
    stratiterate(strat, state) do layer1, layer2, state1, state2
        CryoGrid.interact!(layer1, layer2, state1, state2)
    end
end

@inline function CryoGrid.computefluxes!(strat::Stratigraphy, state)
    fastiterate(layers(strat)) do named_layer
        CryoGrid.computefluxes!(named_layer.val, getproperty(state, layername(named_layer)))
    end
end

@inline function CryoGrid.timestep(strat::Stratigraphy, state)
    max_dts = fastmap(layers(strat)) do named_layer
        CryoGrid.timestep(named_layer.val, getproperty(state, layername(named_layer)))
    end
    return minimum(max_dts)
end

# collecting/grouping components
CryoGrid.events(strat::Stratigraphy) = map(named_layer -> _addlayerfield(CryoGrid.events(named_layer.val), nameof(named_layer)), NamedTuple(strat))

CryoGrid.variables(::Union{FixedVolume,DiagnosticVolume}) = (
    Diagnostic(:Δz, Scalar, u"m", domain=0..Inf),
    # technically the domain for z should be double bounded by the surrounding layers...
    # unfortunately there is no way to represent that here, so we just have to ignore it
    Diagnostic(:z, Scalar, u"m"),
)
CryoGrid.variables(::PrognosticVolume) = (
    Prognostic(:Δz, Scalar, u"m", domain=0..Inf),
    Diagnostic(:z, Scalar, u"m"),
)
# collects all variables in the stratgriphy, returning a NamedTuple of variable sets.
function CryoGrid.variables(strat::Stratigraphy)
    strat_nt = NamedTuple(strat)
    layervars = map(CryoGrid.variables, strat_nt)
    return map(layervars, strat_nt) do vars, named_layer
        (
            vars...,
            CryoGrid.variables(CryoGrid.Volume(named_layer.val))...,
        )
    end
end
# collects and validates all variables in a given layer
function CryoGrid.variables(@nospecialize(named_layer::NamedLayer))
    layer = named_layer.val
    declared_vars = variables(layer)
    nested_vars = Flatten.flatten(layer, Flatten.flattenable, Var)
    all_vars = vcat(collect(declared_vars), collect(nested_vars))
    # check for (permissible) duplicates between variables, excluding parameters
    groups = Utils.groupby(var -> varname(var), all_vars)
    for (id,vargroup) in filter(g -> length(g.second) > 1, groups)
        # if any duplicate variable deifnitions do not match, raise an error
        @assert all(vargroup[i] == vargroup[i-1] for i in 2:length(vargroup)) "Found one or more conflicting definitions of $id in $vargroup"
    end
    diag_vars = filter(isdiagnostic, all_vars)
    prog_vars = filter(isprognostic, all_vars)
    alg_vars = filter(isalgebraic, all_vars)
    # check for duplicated algebraic/prognostic vars by taking intersection of both sets of vars
    prog_alg_duplicated = prog_vars ∩ alg_vars
    @assert isempty(prog_alg_duplicated) "Variables $(prog_alg_duplicated) cannot be both prognostic and algebraic."
    # check for re-definition of diagnostic variables as prognostic;
    # we take the union of both algebraic and prognostic variables and then check for matching diagnostic variables
    prog_alg = prog_vars ∪ alg_vars
    diag_prog = filter(v -> v ∈ prog_alg, diag_vars)
    # check for conflicting definitions of differential vars
    diff_varnames = map(v -> varname(DVar(v)), prog_alg)
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
