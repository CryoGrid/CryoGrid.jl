const RESERVED_LAYER_NAMES = (:top, :bottom, :strat, :init, :event)

# NamedLayer type (alias of Named for Layer types)
const NamedLayer{name,TLayer} = Named{name,TLayer} where {name,TLayer<:Layer}

"""
    Stratigraphy{N,TLayers<:NamedTuple,TBoundaries}

Defines a 1-dimensional stratigraphy by connecting a top and bottom layer to one or more subsurface layers.
"""
struct Stratigraphy{N,TLayers<:NamedTuple,TBoundaries}
    boundaries::TBoundaries
    layers::TLayers
    Stratigraphy(boundaries::NTuple{N,Any}, layers::NamedTuple) where {N} = new{N,typeof(layers),typeof(boundaries)}(boundaries, layers)
    Stratigraphy(
        top::Pair{<:Number,<:Top},
        sub::Pair{<:Number},
        bot::Pair{<:Number,<:Bottom}
    ) = Stratigraphy(top,(sub,),bot)
    Stratigraphy(
        top::Pair{<:Number,<:Top},
        sub::AbstractVector{<:Pair{<:Number}},
        bot::Pair{<:Number,<:Bottom}
    ) = Stratigraphy(top, Tuple(sub), bot)
    Stratigraphy(
        top::Pair{<:Number,<:Top},
        sub::Numerics.Profile{N,<:Tuple{Vararg{Number}},<:Tuple{Vararg{SubSurface}}},
        bot::Pair{<:Number,<:Bottom}
    ) where N = Stratigraphy(top, Tuple(sub), bot)
    function Stratigraphy(
        # use @nospecialize to (hopefully) reduce compilation overhead
        @nospecialize(top::Pair{<:Number,<:Top}),
        @nospecialize(sub::Tuple{Vararg{Pair{<:Number}}}),
        @nospecialize(bot::Pair{<:Number,<:Bottom})
    )
        updateparam(p::Param) = Param(merge(parent(p), (layer=:strat,)))
        updateparam(x) = x
        # check subsurface layers
        @assert length(sub) > 0 "At least one subsurface layer must be specified"
        length(sub) > 18 && @warn "Stratigraphies with more than 20 layers will result in very long compile times and are not recommended. Consider creating heterogeneous layer types."
        top = top[1] => Named(:top => top[2])
        bot = bot[1] => Named(:bottom => bot[2])
        sub = _withnames(sub)
        sub_names = map(nameof, map(last, sub))
        @assert all(map(∉(RESERVED_LAYER_NAMES), sub_names)) "Subsurface layer names may not be one of: $RESERVED_LAYER_NAMES"
        boundaries = map(updateparam ∘ first, (top, sub..., bot))
        # in case one or more boundaries has parameters/units, strip to get numerical value
        boundaryvalues = map(pstrip, boundaries)
        # check boundary depths
        @assert issorted(boundaryvalues) "Stratigraphy boundary locations must be in strictly increasing order."
        # wrap layers
        named_layers = tuple(last(top), map(last, sub)..., last(bot))
        layers = NamedTuple(named_layers)
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
layernames(strat::Stratigraphy) = keys(layers(strat))
layertypes(::Type{<:Stratigraphy{N,NamedTuple{names,TLayers}}}) where {N,TLayers,names} = Tuple(TLayers.parameters)
namedlayers(strat::Stratigraphy) = map(Named, layernames(strat), values(layers(strat)))

Base.keys(strat::Stratigraphy) = layernames(strat)
Base.values(strat::Stratigraphy) = values(layers(strat))
Base.propertynames(strat::Stratigraphy) = Base.keys(strat)
Base.getproperty(strat::Stratigraphy, name::Symbol) = getproperty(layers(strat), name)
Base.getindex(strat::Stratigraphy, name::Symbol) = getproperty(strat, name)
# Array and iteration overrides
Base.size(strat::Stratigraphy) = size(layers(strat))
Base.length(strat::Stratigraphy) = length(layers(strat))
Base.getindex(strat::Stratigraphy, i::Int) = layers(strat)[i]
Base.firstindex(strat::Stratigraphy) = 1
Base.lastindex(strat::Stratigraphy) = length(strat)
Base.iterate(strat::Stratigraphy) = iterate(layers(strat))
Base.iterate(strat::Stratigraphy, i) = iterate(layers(strat), i)
Base.NamedTuple(strat::Stratigraphy) = layers(strat)
function Base.show(io::IO, ::MIME"text/plain", strat::Stratigraphy)
    print(io, "Stratigraphy:\n")
    for (k,v,b) in zip(keys(strat), values(strat), boundaries(strat))
        print(repeat(" ", Base.indent_width), "$b: $k :: $(nameof(typeof(v)))\n")
    end
end

function Numerics.makegrid(strat::Stratigraphy, strategy::DiscretizationStrategy)
    strat_grid = nothing
    for (bounds, named_layer) in zip(boundarypairs(strat)[2:end-1], namedlayers(strat)[2:end-1])
        @assert bounds[2] - bounds[1] > zero(bounds[1]) "Subsurface layers must have thickness greater than zero in the stratigraphy. The initial thickness can be later updated in `initialcondition!`."
        layer = named_layer
        layer_grid = Numerics.makegrid(layer.val, strategy, bounds)
        if !isnothing(strat_grid)
            # check that grid edges line up at layer boundary
            @assert strat_grid[end] == layer_grid[1] "Upper boundary of layer $(nameof(named_layer)) does not match the previous layer."
            # concatenate grids, omitting first value of layer_grid to avoid duplicating the shared edge
            strat_grid = Grid(vcat(strat_grid, layer_grid[2:end]))
        else
            strat_grid = layer_grid
        end
    end
    return Grid(ustrip.(strat_grid))
end

function CryoGrid.initialcondition!(strat::Stratigraphy, state, inits::CryoGrid.Initializer...)
    # initialcondition! is only called once so we don't need to worry about performance;
    # we can just loop over everything.

    varinits = filter(init -> isa(init, VarInitializer), inits)
    otherinits = filter(init -> !isa(init, VarInitializer), inits)

    # first invoke var initializers
    _initializers!(strat, state, varinits)

    # then invoke non-specific layer inits
    for i in 1:length(strat)-1
        layerᵢ = strat[i]
        layerᵢ₊₁ = strat[i+1]
        stateᵢ = getproperty(state, layernames(strat)[i])
        stateᵢ₊₁ = getproperty(state, layernames(strat)[i+1])
        if i == 1
            CryoGrid.initialcondition!(layerᵢ, stateᵢ)
        end
        CryoGrid.initialcondition!(layerᵢ₊₁, stateᵢ₊₁)
        CryoGrid.initialcondition!(layerᵢ, layerᵢ₊₁, stateᵢ, stateᵢ₊₁)
    end

    # and finally non-var initializers
    _initializers!(strat, state, otherinits)
end

function CryoGrid.resetfluxes!(strat::Stratigraphy, state)
    fastiterate(namedlayers(strat)) do named_layer
        CryoGrid.resetfluxes!(named_layer.val, getproperty(state, nameof(named_layer)))
    end
end

function CryoGrid.computediagnostic!(strat::Stratigraphy, state)
    fastiterate(namedlayers(strat)) do named_layer
        CryoGrid.computediagnostic!(named_layer.val, getproperty(state, nameof(named_layer)))
    end
end

function CryoGrid.diagnosticstep!(strat::Stratigraphy, state)
    prognostic_state_updated = fastmap(namedlayers(strat)) do named_layer
        CryoGrid.diagnosticstep!(named_layer.val, getproperty(state, nameof(named_layer)))
    end
    return any(prognostic_state_updated)
end

"""
    interact!(strat::Stratigraphy, state)

Special implementation of `interact!` that iterates over each pair of layers in the stratigraphy which are adjacent and
"active" based on the current `state`. `state` must have properties defined corresponding to the name of each layer such
that `getproperty` would return the appropriate state object for the i'th layer in the stratigraphy.
"""
@generated function CryoGrid.interact!(strat::Stratigraphy{N}, state) where {N}
    expr = Expr(:block)
    # build expressions for checking whether each layer is active
    is_active_exprs = map(i -> :(CryoGrid.isactive(strat[$i], getproperty(state, layernames(strat)[$i]))), tuple(1:N...))
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
                    # always try to invoke interact! for adjacent layers
                    quote
                        if caninteract(strat[$i], strat[$j], getproperty(state, names[$i]), getproperty(state, names[$j]))
                            interact!(strat[$i], strat[$j], getproperty(state, names[$i]), getproperty(state, names[$j]))
                        end
                    end
                else
                    # for layers that are non-adjacent, we only invoke interact! if all layers between them are currently inactive.
                    inactive_imbetween_layers = map(i+1:j-1) do k
                        quote
                            !is_active[$k]
                        end
                    end
                    quote
                        # if layers i and j are both active, and all imbetween layers are inactive, invoke f!;
                        # note the splat syntax here: $(inactive_imbetween_layers...) simply expands the tuple of
                        # expressions 'inactive_imbetween_layers' as arguments to `tuple`.
                        can_interact = caninteract(strat[$i], strat[$j], getproperty(state, names[$i]), getproperty(state, names[$j]))
                        if is_active[$i] && is_active[$j] && all(tuple($(inactive_imbetween_layers...))) && can_interact
                            interact!(strat[$i], strat[$j], getproperty(state, names[$i]), getproperty(state, names[$j]))
                        end
                    end
                end
            push!(expr.args, expr_ij)
        end
    end
    push!(expr.args, :(return nothing))
    return expr
end

function CryoGrid.computefluxes!(strat::Stratigraphy, state)
    fastiterate(namedlayers(strat)) do named_layer
        CryoGrid.computefluxes!(named_layer.val, getproperty(state, nameof(named_layer)))
    end
end

function CryoGrid.timestep(strat::Stratigraphy, state)
    @inline function timestep(named_layer::NamedLayer)
        return CryoGrid.timestep(named_layer.val, getproperty(state, nameof(named_layer)))
    end
    max_dts = fastmap(timestep, namedlayers(strat))
    return minimum(max_dts)
end

# collecting/grouping components
CryoGrid.events(strat::Stratigraphy) = (; map(named_layer -> nameof(named_layer) => _addlayerfield(CryoGrid.events(named_layer.val), nameof(named_layer)), namedlayers(strat))...)

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
    layervars = map(_collectvars, layers(strat))
    return map(layervars, layers(strat)) do vars, layer
        (
            vars...,
            CryoGrid.variables(CryoGrid.Volume(layer))...,
        )
    end
end

function _initializers!(strat::Stratigraphy, state, inits)
    for i in 1:length(strat)
        layerᵢ = strat[i]
        stateᵢ = getproperty(state, layernames(strat)[i])
        if i < length(strat)
            layerᵢ₊₁ = strat[i+1]
            stateᵢ₊₁ = getproperty(state, layernames(strat)[i+1])
        end
        # loop over initializers
        for init in sort(collect(tuplejoin(CryoGrid.initializers(layerᵢ), inits)))
            isvalidᵢ = !isa(init, VarInitializer) || hasproperty(stateᵢ, varname(init))
            if isvalidᵢ
                CryoGrid.initialcondition!(init, layerᵢ, stateᵢ)
            end

            if i < length(strat)
                isvalidᵢ₊₁ = !isa(init, VarInitializer) || hasproperty(stateᵢ₊₁, varname(init))
                if isvalidᵢ && isvalidᵢ₊₁
                    CryoGrid.initialcondition!(init, layerᵢ, layerᵢ₊₁, stateᵢ, stateᵢ₊₁)
                end
            end
        end
    end
end

# collects and validates all variables in a given layer
function _collectvars(@nospecialize(layer::Layer))
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

function _withnames(pairs::Tuple{Vararg{Pair}})
    # case 1: user specified name as :name => layer or Named
    withname(l::Pair{Symbol,<:Layer}) = Named(l)
    withname(l::NamedLayer) = l
    # case 2: no name specified, autogenerate based on type name
    withname(l::Layer) = Named(Symbol(lowercase(string(nameof(typeof(l))))), l)
    # extract depths and layers
    depths = map(first, pairs)
    layers = map(last, pairs)
    # get named layers
    named_layers = map(withname, layers)
    # group by name
    grouped_layers = Utils.groupby(nameof, named_layers)
    # iterate over each name group
    for name in keys(grouped_layers)
        if length(grouped_layers[name]) > 1
            # if there is more than one layer with this name, append an integer to deduplicate them.
            grouped_layers[name] = map(enumerate(grouped_layers[name])) do (i, named_layer)
                # rename layer with integer id to deduplicate
                Named(Symbol(nameof(named_layer), i), named_layer.val)
            end
        end
    end
    return map(depths, named_layers) do depth, named_layer
        depth => popfirst!(grouped_layers[nameof(named_layer)])
    end
end
