mutable struct StateHistory
    vals::Union{Missing,<:Any}
    StateHistory() = new(missing)
end

"""
    AbstractLandModel{iip}

Base type for 1D land models. `iip` is a value of enum `InPlaceMode` that indicates
whether the model operates on state variables in-place (overwriting arrays) or
out-of-place (copying arrays).
"""
abstract type AbstractLandModel{iip} end
"""
    (model::AbstractLandModel{inplace})(du,u,p,t)
    (model::AbstractLandModel{ooplace})(u,p,t)

Invokes the corresponding `step` function to compute the time derivative du/dt.
"""
(model::AbstractLandModel{inplace})(du,u,p,t) = step!(model,du,u,p,t)
(model::AbstractLandModel{ooplace})(u,p,t) = step(model,u,p,t)

"""
    step!(::T,du,u,p,t) where {T<:AbstractLandModel}

In-place step function for model `T`. Computes du/dt and stores the result in `du`.
"""
step!(::T,du,u,p,t) where {T<:AbstractLandModel} = error("no implementation of in-place step! for $T")
"""
    step(::T,u,p,t) where {T<:AbstractLandModel}

Out-of-place step function for model `T`. Computes and returns du/dt as vector with same size as `u`.
"""
step(::T,u,p,t) where {T<:AbstractLandModel} = error("no implementation of out-of-place step for $T")

"""
    LandModel{TStrat,TGrid,TStates,iip,obsv} <: AbstractLandModel{iip}

Defines the full specification of a CryoGrid land model; i.e. stratigraphy, grid, and state variables.
"""
struct LandModel{TStrat,TGrid,TStates,iip,obsv} <: AbstractLandModel{iip}
    strat::TStrat # stratigraphy
    grid::TGrid # grid
    state::TStates # state variables
    hist::StateHistory # mutable "history" type for state tracking
    function LandModel(
        strat::TStrat,
        grid::TGrid,
        state::TStates,
        hist::StateHistory=StateHistory(),
        iip::InPlaceMode=inplace,
        observe::Vector{Symbol}=Symbol[]) where
        {TStrat<:Stratigraphy,TGrid<:Grid{Edges},TStates<:VarStates}
        new{TStrat,TGrid,TStates,iip,tuple(observe...)}(strat,grid,state,hist)
    end
end
ConstructionBase.constructorof(::Type{LandModel{TStrat,TGrid,TStates,iip,obsv}}) where {TStrat,TGrid,TStates,iip,obsv} =
    (strat, grid, state, hist) -> LandModel(strat,grid,state,hist,iip,length(obsv) > 0 ? collect(obsv) : Symbol[])

Base.show(io::IO, ::MIME"text/plain", model::LandModel{TStrat,TGrid,TStates,iip,obsv}) where {TStrat,TGrid,TStates,iip,obsv} = print(io, "LandModel ($iip) with layers $(map(componentname, components(model.strat))), observables=$obsv, $TGrid, $TStrat")

"""
Constructs a `LandModel` from the given stratigraphy and grid. `arrayproto` keyword arg should be an array instance
(of any arbitrary length, including zero, contents are ignored) that will determine the array type used for all state vectors.
"""
function LandModel(
    @nospecialize(strat::Stratigraphy),
    @nospecialize(grid::Grid{Edges,<:Numerics.Geometry,<:DistQuantity});
    arrayproto::Type{A}=Vector,
    iip::InPlaceMode=inplace,
    observe::Vector{Symbol}=Symbol[],
    chunksize=nothing,
) where {A<:AbstractArray}
    vars = OrderedDict()
    components = OrderedDict()
    for (i,comp) in enumerate(strat)
        name = componentname(comp)
        # build layer
        vars[name] = _collectvars(comp)
        params = ModelParameters.params(comp)
        if length(params) > 0
            # create sub-model and add layer name to all parameters
            submodel = Model(comp)
            submodel[:layer] = repeat([name], length(params))
            components[name] = parent(submodel)
        else
            components[name] = comp
        end
    end
    # rebuild stratigraphy with updated parameters
    strat = Stratigraphy(boundaries(strat), Tuple(values(components)))
    para = params(strat)
    # construct state variables
    componentnames = [componentname(node) for node in strat]
    ntvars = NamedTuple{Tuple(componentnames)}(Tuple(values(vars)))
    npvars = (length(filter(isprognostic, layer)) + length(filter(isalgebraic, layer)) for layer in ntvars) |> sum
    ndvars = (length(filter(isdiagnostic, layer)) for layer in ntvars) |> sum
    @assert (npvars + ndvars) > 0 "No variable definitions found. Did you add a method definition for CryoGrid.variables(::L,::P) where {L<:Layer,P<:Process}?"
    @assert npvars > 0 "At least one prognostic variable must be specified."
    chunksize = isnothing(chunksize) ? length(para) : chunksize
    states = VarStates(ntvars, Grid(dustrip(grid), grid.geometry), chunksize, arrayproto)
    LandModel(strat,grid,states,StateHistory(),iip,observe)
end
LandModel(strat::Stratigraphy, grid::Grid{Cells}; kwargs...) = LandModel(strat, edges(grid); kwargs...)
LandModel(strat::Stratigraphy, grid::Grid{Edges,<:Numerics.Geometry,T}; kwargs...) where {T} = error("grid must have values with units of length, e.g. try using `Grid((x)u\"m\")` where `x` are your grid points.")
# mark only stratigraphy field as flattenable
Flatten.flattenable(::Type{<:LandModel}, ::Type{Val{:strat}}) = true
Flatten.flattenable(::Type{<:LandModel}, ::Type{Val{name}}) where name = false

"""
Generated step function (i.e. du/dt) for any arbitrary LandModel. Specialized code is generated and compiled
on the fly via the @generated macro to ensure type stability. The generated code updates each layer in the stratigraphy
in sequence, i.e for each layer 1 < i < N:

diagnosticstep!(layer i, ...)
interact!(layer i-1, ...)
prognosticstep!(layer i, ...)

Note for developers: All sections of code wrapped in quote..end blocks are generated. Code outside of quote blocks
is only executed during compilation and will not appear in the compiled version.
"""
@generated function step!(model::LandModel{TStrat,TGrid,TStates,inplace,obsv}, _du,_u,_p,t) where {TStrat,TGrid,TStates,obsv}
    nodetyps = componenttypes(TStrat)
    N = length(nodetyps)
    expr = Expr(:block)
    # Declare variables
    @>> quote
    p = updateparams!(_p, model, _du, _u, t)
    strat = Flatten.reconstruct(model.strat, p, ModelParameters.SELECT, ModelParameters.IGNORE)
    _du .= zero(eltype(_du))
    du = ComponentArray(_du, getaxes(model.state.uproto))
    u = ComponentArray(_u, getaxes(model.state.uproto))
    state = LandModelState(model.state, boundaries(strat), u, du, t, Val{inplace}())
    end push!(expr.args)
    # Initialize variables for all layers
    for i in 1:N
        n = componentname(nodetyps[i])
        nstate = Symbol(n,:state)
        nlayer = Symbol(n,:layer)
        nprocess = Symbol(n,:process)
        @>> quote
        $nstate = state.$n
        $nlayer = components(strat)[$i].layer
        $nprocess = components(strat)[$i].process
        end push!(expr.args)
    end
    # Diagnostic step
    for i in 1:N
        n = componentname(nodetyps[i])
        nstate = Symbol(n,:state)
        nlayer = Symbol(n,:layer)
        nprocess = Symbol(n,:process)
        @>> quote
        diagnosticstep!($nlayer,$nprocess,$nstate)
        end push!(expr.args)
    end
    # Interact
    for i in 1:N-1
        n1,n2 = componentname(nodetyps[i]), componentname(nodetyps[i+1])
        n1state, n2state = Symbol(n1,:state), Symbol(n2,:state)
        n1layer, n2layer = Symbol(n1,:layer), Symbol(n2,:layer)
        n1process, n2process = Symbol(n1,:process), Symbol(n2,:process)
        @>> quote
        interact!($n1layer,$n1process,$n2layer,$n2process,$n1state,$n2state)
        end push!(expr.args)
    end
    # Prognostic step
    for i in 1:N
        n = componentname(nodetyps[i])
        nstate = Symbol(n,:state)
        nlayer = Symbol(n,:layer)
        nprocess = Symbol(n,:process)
        @>> quote
        prognosticstep!($nlayer,$nprocess,$nstate)
        end push!(expr.args)
    end
    # Observables
    for i in 1:N
        n = componentname(nodetyps[i])
        nstate = Symbol(n,:state)
        nlayer = Symbol(n,:layer)
        nprocess = Symbol(n,:process)
        for name in obsv
            nameval = Val{name}()
            @>> quote
                observe($nameval,$nlayer,$nprocess,$nstate)
            end push!(expr.args)
        end
    end
    # make sure compiled method returns no value
    @>> :(return nothing) push!(expr.args)
    # emit generated expression block
    return expr
end
"""
    initialcondition!(model::LandModel, tspan::NTuple{2,Float64}, p::AbstractVector, initializers::VarInit...)
    initialcondition!(model::LandModel, tspan::NTuple{2,DateTime}, p::AbstractVector, initializers::VarInit...)

Calls `initialcondition!` on all layers/processes and returns the fully constructed u0 and du0 states.
"""
initialcondition!(model::LandModel, tspan::NTuple{2,DateTime}, p::AbstractVector, args...) = initialcondition!(model, convert_tspan(tspan), p, args...)
@generated function initialcondition!(model::LandModel{TStrat,TGrid,TStates,iip,obsv}, tspan::NTuple{2,Float64}, p::AbstractVector, initializers::Numerics.VarInit...) where {TStrat,TGrid,TStates,iip,obsv}
    nodetyps = componenttypes(TStrat)
    N = length(nodetyps)
    expr = Expr(:block)
    # Declare variables
    @>> quote
    du = zero(similar(model.state.uproto, eltype(p)))
    u = zero(similar(model.state.uproto, eltype(p)))
    strat = Flatten.reconstruct(model.strat, p, ModelParameters.SELECT, ModelParameters.IGNORE)
    state = LandModelState(model.state, boundaries(strat), u, du, tspan[1], Val{iip}())
    end push!(expr.args)
    # Call initializers
    for i in 1:N
        for j in 1:length(initializers)
            @>> quote
            let layerstate = state[$i],
                init = initializers[$j];
                if haskey(layerstate.states, varname(init))
                    initvar!(layerstate, strat, init)
                end
            end
            end push!(expr.args)
        end
    end
    # Iterate over layers
    for i in 1:N-1
        n1,n2 = componentname(nodetyps[i]), componentname(nodetyps[i+1])
        # create variable names
        n1,n2 = componentname(nodetyps[i]), componentname(nodetyps[i+1])
        n1state, n2state = Symbol(n1,:state), Symbol(n2,:state)
        n1layer, n2layer = Symbol(n1,:layer), Symbol(n2,:layer)
        n1process, n2process = Symbol(n1,:process), Symbol(n2,:process)
        # generated code for layer updates
        @>> quote
        $n1state = state.$n1
        $n2state = state.$n2
        $n1layer = components(strat)[$i].layer
        $n1process = components(strat)[$i].process
        $n2layer = components(strat)[$(i+1)].layer
        $n2process = components(strat)[$(i+1)].process
        end push!(expr.args)
        if i == 1
            # only invoke initialcondition! for layer i on first iteration to avoid duplicated calls
            @>> quote
            initialcondition!($n1layer,$n1state)
            initialcondition!($n1layer,$n1process,$n1state)
            end push!(expr.args)
        end
        # invoke initialcondition! for each layer, then for both (similar to interact!)
        @>> quote
        initialcondition!($n2layer,$n2state)
        initialcondition!($n2layer,$n2process,$n2state)
        initialcondition!($n1layer,$n1process,$n2layer,$n2process,$n1state,$n2state)
        end push!(expr.args)
    end
    @>> quote
    return u
    end push!(expr.args)
    return expr
end
"""
    initvar!(state::LayerState, ::Stratigraphy, init::VarInit{varname}) where {varname}
    initvar!(state::LayerState, ::Stratigraphy, init::InterpInit{varname})

Calls the initializer for state variable `varname`.
"""
initvar!(state::LayerState, ::Stratigraphy, init::Numerics.VarInit{varname}) where {varname} = init!(state[varname], init)
initvar!(state::LayerState, ::Stratigraphy, init::Numerics.InterpInit{varname}) where {varname} = init!(state[varname], init, state.grids[varname])
"""
    getvar(name::Symbol, model::LandModel, u)
    getvar(::Val{name}, model::LandModel, u)

Retrieves the (diagnostic or prognostic) grid variable from `model` given prognostic state `u`.
If `name` is not a variable in the model, or if it is not a grid variable, `nothing` is returned.
"""
Numerics.getvar(name::Symbol, model::LandModel, u) = getvar(Val{name}(), model, u)
Numerics.getvar(::Val{name}, model::LandModel, u) where name = getvar(Val{name}(), model.state, withaxes(u, model))
"""
    getstate(layername::Symbol, model::LandModel, u, du, t)
    getstate(::Val{layername}, model::LandModel{TStrat,TGrid,<:VarStates{layernames},iip}, _u, _du, t)

Constructs a `LayerState` representing the full state of `layername` given `model`, state vectors `u` and `du`, and the
time step `t`.
"""
getstate(layername::Symbol, model::LandModel, u, du, t) = getstate(Val{layername}(), model, u, du, t)
function getstate(::Val{layername}, model::LandModel{TStrat,TGrid,<:VarStates{layernames},iip}, _u, _du, t) where {layername,TStrat,TGrid,iip,layernames}
    du = ComponentArray(_du, getaxes(model.state.uproto))
    u = ComponentArray(_u, getaxes(model.state.uproto))
    i = 1
    for j in 1:length(model.strat)
        if layernames[j] == layername
            i = j
            break
        end
    end
    z = boundaryintervals(map(ustrip, stripparams(boundaries(model.strat))), ustrip(model.grid[end]))[i]
    return LayerState(model.state, z, u, du, t, Val{layername}(), Val{iip}())
end
"""
    variables(model::LandModel)

Returns a tuple of all variables defined in the model.
"""
variables(model::LandModel) = Tuple(unique(Flatten.flatten(model.state.vars, Flatten.flattenable, Var)))
"""
    withaxes(u::AbstractArray, ::LandModel)

Constructs a `ComponentArray` with labeled axes from the given state vector `u`. Assumes `u` to be of the same type/shape
as `setup.uproto`.
"""
withaxes(u::AbstractArray, model::LandModel) = ComponentArray(u, getaxes(model.state.uproto))
withaxes(u::ComponentArray, ::LandModel) = u
"""
Gets the 
"""
function getstate(model::LandModel{TStrat,TGrid,TStates,iip}, _u, _du, t) where {TStrat,TGrid,TStates,iip}
    du = ComponentArray(_du, getaxes(model.state.uproto))
    u = ComponentArray(_u, getaxes(model.state.uproto))
    return LandModelState(model.strat, model.state, u, du, t, Val{iip}())
end
"""
Collects and validates all declared variables (`Var`s) for the given strat component.
"""
function _collectvars(@nospecialize(comp::StratComponent))
    layer, process = comp.layer, comp.process
    layer_vars = variables(layer)
    @assert all([isdiagnostic(var) for var in layer_vars]) "Layer variables must be diagnostic."
    process_vars = variables(layer, process)
    all_vars = tuple(layer_vars...,process_vars...)
    @debug "Building layer $(componentname(comp)) with $(length(all_vars)) variables: $(all_vars)"
    # check for (permissible) duplicates between variables, excluding parameters
    groups = groupby(var -> varname(var), all_vars)
    for (name,gvars) in filter(g -> length(g.second) > 1, groups)
        # if any duplicate variable deifnitions do not match, raise an error
        @assert reduce(==, gvars) "Found one or more conflicting definitions of $name in $gvars"
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
    if !isempty(diag_prog)
        @warn "Variables $(Tuple(map(varname,diag_prog))) declared as both prognostic/algebraic and diagnostic. In-place modifications outside of callbacks may degrade integration accuracy."
    end
    # check for conflicting definitions of differential vars
    diff_varnames = map(v -> Symbol(:d, varname(v)), prog_alg)
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
