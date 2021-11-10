using DiffEqBase
using DiffEqCallbacks
using OrdinaryDiffEq

@reexport using DiffEqBase: solve, init, ODEProblem, SciMLBase

"""
Specialized problem type for CryoGrid `ODEProblem`s.
"""
struct CryoGridODEProblem end

"""
    CryoGridProblem(setup::LandModel, tspan::NTuple{2,Float64}, p=nothing;kwargs...)

CryoGrid specialized constructor for ODEProblem that automatically generates the initial
condition and necessary callbacks.
"""
function CryoGridProblem(setup::LandModel, tspan::NTuple{2,Float64}, p=nothing; saveat=3600.0, save_everystep=false, callback=nothing, save_adtypes=false, kwargs...)
    # workaround for bug in DiffEqCallbacks; see https://github.com/SciML/DifferentialEquations.jl/issues/326
    # we have to manually expand single-number `saveat` (i.e. time interval for saving) to a step-range.
    expandtstep(tstep::Number) = tspan[1]:tstep:tspan[end]
    expandtstep(tstep::AbstractVector) = tstep
    values(u, ::Val{false}) = reinterpret(Float64, u) # the actual values resulting from reinterpret here will be garbage, but we only need the type for `getstates`
    values(u, ::Val{true}) = true
    model = Model(setup)
    p = isnothing(p) ? dustrip.(collect(model[:val])) : p
	# compute initial condition
	u0, du0 = Land.init!(setup, tspan, p)
    # set up saving callback
    stateproto = getstates(setup, values(du0,Val{save_adtypes}()), values(u0,Val{save_adtypes}()), tspan[1], Val{:diagnostic}())
    savevals = SavedValues(Float64, typeof(stateproto))
    savefunc = (u,t,integrator) -> deepcopy(getstates(setup, values(get_du(integrator),Val{save_adtypes}()), values(u,Val{save_adtypes}()), t, Val{:diagnostic}()))
    savingcallback = SavingCallback(savefunc, savevals; saveat=expandtstep(saveat), save_everystep=save_everystep)
    callbacks = isnothing(callback) ? savingcallback : CallbackSet(savingcallback, callback)
    # note that this implicitly discards any existing saved values in the model setup's state history
    setup.hist.vals = savevals
    # set up default mass matrix
    M_diag = similar(setup.uproto)
    for layer in keys(setup.meta)
        progvars = setup.meta[layer][:progvars]
        algvars = setup.meta[layer][:algvars]
        M_diag_layer = @view M_diag[layer]
        for var in progvars
            M_diag_layer[varname(var)] .= one(eltype(M_diag))
        end
        for var in algvars
            M_diag_layer[varname(var)] .= zero(eltype(M_diag))
        end
    end
    # if no algebraic variables are present, use identity matrix
    num_algebraic = length(M_diag) - sum(M_diag)
    M = num_algebraic > 0 ? Diagonal(M_diag) : I
	func = odefunction(setup, M, u0, p, tspan; kwargs...)
	ODEProblem(func, u0, tspan, p, CryoGridODEProblem(); callback=callbacks, kwargs...)
end
# this version converts tspan from DateTime to float
"""
    CryoGridProblem(setup::LandModel, tspan::NTuple{2,DateTime}, args...;kwargs...)
"""
CryoGridProblem(setup::LandModel, tspan::NTuple{2,DateTime}, args...;kwargs...) = CryoGridProblem(setup,convert_tspan(tspan),args...;kwargs...)

export CryoGridProblem

"""
    odefunction(setup::LandModel, M, u0, p, tspan; kwargs...)

Constructs a SciML `ODEFunction` given the model setup, mass matrix M, initial state u0, parameters p, and tspan.
Can (and should) be overridden by users to provide customized ODEFunction configurations for specific problem setups, e.g:
```
model = LandModel(strat,grid)
function CryoGrid.Setup.odefunction(::DefaultJac, setup::typeof(model), M, u0, p, tspan)
    ...
    # make sure to return an instance of ODEFunction
end
...
prob = CryoGridProblem(model, tspan, p)
```

`JacobianStyle` can also be extended to create custom traits which can then be applied to compatible `LandModel`s.
"""
odefunction(setup::TSetup, M, u0, p, tspan; kwargs...) where {TSetup<:LandModel} = odefunction(JacobianStyle(TSetup), setup, M, u0, p, tspan; kwargs...)
odefunction(::DefaultJac, setup::TSetup, M, u0, p, tspan; kwargs...) where {TSetup<:LandModel} = ODEFunction(setup, mass_matrix=M; kwargs...)
function odefunction(::TridiagJac, setup::LandModel, M, u0, p, tspan; kwargs...)
    if :jac_prototype in keys(kwargs)
        @warn "using user specified jac_prorotype instead of tridiagonal"
        ODEFunction(setup, mass_matrix=M; kwargs...)
    else
        N = length(u0)
        J = Tridiagonal(
                similar(u0, eltype(p), N-1) |> Vector,
                similar(u0, eltype(p), N) |> Vector,
                similar(u0, eltype(p), N-1) |> Vector
        )
        ODEFunction(setup, mass_matrix=M, jac_prototype=J, kwargs...)
    end
end

"""
    getstate(layername::Symbol, integrator::SciMLBase.DEIntegrator)

Builds the state named tuple for `layername` given an initialized integrator.
"""
getstate(layername::Symbol, integrator::SciMLBase.DEIntegrator) = getstate(Val{layername}(), integrator)
@generated function getstate(::Val{layername}, integrator::SciMLBase.DEIntegrator) where {layername}
    # a bit hacky and may break in the future... but this is the hardcoded position of the LandModel type in DEIntegrator
    TStrat = integrator.parameters[13].parameters[2].parameters[1]
    names = map(componentname, componenttypes(TStrat))
    i = findfirst(n -> n == layername, names)
    quote
        let setup = integrator.f.f;
            Land._buildstate(
                setup.cache[$(QuoteNode(layername))],
                setup.meta[$(QuoteNode(layername))],
                withaxes(integrator.u,setup).$layername,
                withaxes(get_du(integrator),setup).$layername,
                integrator.t,
                boundaries(setup.strat)[$i]
            )
        end
    end
end
"""
    getvar(var::Symbol, integrator::SciMLBase.DEIntegrator)
"""
getvar(var::Symbol, integrator::SciMLBase.DEIntegrator) = Land.getvar(Val{var}(), integrator.f.f, integrator.u)
