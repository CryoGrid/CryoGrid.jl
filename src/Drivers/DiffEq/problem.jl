"""
Specialized problem type for CryoGrid `ODEProblem`s.
"""
struct CryoGridODEProblem end

"""
    CryoGridProblem(setup::Tile, tspan::NTuple{2,Float64}, p=nothing;kwargs...)

CryoGrid specialized constructor for ODEProblem that automatically generates the initial
condition and necessary callbacks.
"""
function CryoGridProblem(
    tile::Tile,
    u0::ComponentVector,
    tspan::NTuple{2,Float64},
    p=nothing;
    saveat=3600.0,
    savevars=(),
    save_everystep=false,
    save_start=true,
    save_end=true,
    step_limiter=CryoGrid.timestep,
    safety_factor=1,
    max_step=true,
    callback=nothing,
    isoutofdomain=Strat.domain(tile),
    specialization=SciMLBase.AutoSpecialize,
    function_kwargs=(),
    prob_kwargs...
)
    # workaround for bug in DiffEqCallbacks; see https://github.com/SciML/DifferentialEquations.jl/issues/326
    # we have to manually expand single-number `saveat` (i.e. time interval for saving) to a step-range.
    expandtstep(tstep::Number) = tspan[1]:tstep:tspan[end]
    expandtstep(tstep::AbstractVector) = tstep
    getsavestate(tile::Tile, u, du) = deepcopy(Strat.getvars(tile.state, Strat.withaxes(u, tile), Strat.withaxes(du, tile), savevars...))
    savefunc(u, t, integrator) = getsavestate(Tile(integrator), Strat.withaxes(u, Tile(integrator)), get_du(integrator))
    model_tile = Model(tile)
    if !isnothing(p)
        model_tile[:val] = p
    end
    tile = parent(model_tile)
    # collect parameters
    p = collect(model_tile[:val])
    du0 = zero(u0)
    # remove units
    tile = stripunits(tile)
    # set up saving callback
    stateproto = getsavestate(tile, u0, du0)
    savevals = SavedValues(Float64, typeof(stateproto))
    savingcallback = SavingCallback(savefunc, savevals; saveat=expandtstep(saveat), save_start=save_start, save_end=save_end, save_everystep=save_everystep)
    # add step limiter to default callbacks, if defined
    defaultcallbacks = isnothing(step_limiter) ? (savingcallback,) : (savingcallback, StepsizeLimiter(step_limiter; safety_factor, max_step))
    # build layer callbacks
    layercallbacks = DiffEq.makecallbacks(tile)
    # add user callbacks
    usercallbacks = isnothing(callback) ? () : callback
    callbacks = CallbackSet(defaultcallbacks..., layercallbacks..., usercallbacks...)
    # note that this implicitly discards any existing saved values in the model setup's state history
    tile.hist.vals = savevals
    # set up default mass matrix, M:
    # M⋅∂u∂t = f(u)
    M_diag = similar(tile.state.uproto)
    M_idxmap = ComponentArrays.indexmap(getaxes(M_diag)[1])
    allvars = Flatten.flatten(tile.state.vars, Flatten.flattenable, Var)
    progvars = map(varname, filter(isprognostic, allvars))
    algvars = map(varname, filter(isalgebraic, allvars))
    for name in keys(M_idxmap)
        M_diag_var = @view M_diag[name]
        if name ∈ progvars
            M_diag_var .= one(eltype(M_diag))
        elseif name ∈ algvars
            M_diag_var .= zero(eltype(M_diag))
        end
    end
    # if no algebraic variables are present, use identity matrix
    num_algebraic = length(M_diag) - sum(M_diag)
    M = num_algebraic > 0 ? Diagonal(M_diag) : I
	func = odefunction(tile, u0, p, tspan; mass_matrix=M, specialization, function_kwargs...)
	ODEProblem(func, u0, tspan, p, CryoGridODEProblem(); callback=callbacks, isoutofdomain, prob_kwargs...)
end
"""
    CryoGridProblem(setup::Tile, tspan::NTuple{2,DateTime}, args...;kwargs...)
"""
CryoGridProblem(setup::Tile, u0::ComponentVector, tspan::NTuple{2,DateTime}, args...;kwargs...) = CryoGridProblem(setup,u0,convert_tspan(tspan),args...;kwargs...)
"""
    odefunction(setup::Tile, u0, p, tspan; kwargs...)

Constructs a SciML `ODEFunction` given the model setup, initial state u0, parameters p, and tspan.
Can (and should) be overridden by users to provide customized ODEFunction configurations for specific problem setups, e.g:
```
tile = Tile(strat,grid)
function CryoGrid.Setup.odefunction(::DefaultJac, setup::typeof(tile), u0, p, tspan)
    ...
    # make sure to return an instance of ODEFunction
end
...
prob = CryoGridProblem(tile, tspan, p)
```

`JacobianStyle` can also be extended to create custom traits which can then be applied to compatible `Tile`s.
"""
odefunction(tile::TTile, u0, p, tspan; mass_matrix=I, specialization=SciMLBase.AutoSpecialize, kwargs...) where {TTile<:Tile} = odefunction(JacobianStyle(TTile), tile, u0, p, tspan; mass_matrix, specialization, kwargs...)
odefunction(::DefaultJac, tile::Tile, u0, p, tspan; mass_matrix=I, specialization=SciMLBase.AutoSpecialize, kwargs...) = ODEFunction{true,specialization}(tile; mass_matrix, kwargs...)
function odefunction(::TridiagJac, tile::Tile, u0, p, tspan; mass_matrix=I, specialization=SciMLBase.AutoSpecialize, kwargs...)
    if :jac_prototype in keys(kwargs)
        @warn "using user specified jac_prorotype instead of tridiagonal"
        ODEFunction{true,specialization}(tile; mass_matrix, kwargs...)
    else
        N = length(u0)
        J = Tridiagonal(
                similar(u0, eltype(p), N-1) |> Vector,
                similar(u0, eltype(p), N) |> Vector,
                similar(u0, eltype(p), N-1) |> Vector
        )
        ODEFunction{true,specialization}(tile; jac_prototype=J, mass_matrix, kwargs...)
    end
end
