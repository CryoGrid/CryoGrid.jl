"""
Specialized problem type for CryoGrid `ODEProblem`s.
"""
struct CryoGridODEProblem end

"""
    CryoGridProblem(setup::CryoGridSetup, tspan::NTuple{2,Float64}, p=nothing;kwargs...)

CryoGrid specialized constructor for ODEProblem that automatically generates the initial
condition and necessary callbacks.
"""
function CryoGridProblem(setup::CryoGridSetup, tspan::NTuple{2,Float64}, p=nothing; saveat=3600.0, save_everystep=false, callback=nothing, kwargs...)
    # workaround for bug in DiffEqCallbacks; see https://github.com/SciML/DifferentialEquations.jl/issues/326
    # we have to manually expand single-number `saveat` (i.e. time interval for saving) to a step-range.
    expandtstep(tstep::Number) = tspan[1]:tstep:tspan[end]
    expandtstep(tstep::AbstractVector) = tstep
	p = isnothing(p) ? setup.pproto : p
	# compute initial condition
	u0, du0 = init!(setup, p, tspan)
    # set up saving callback
    stateproto = getstates(setup, du0, u0, p, tspan[1], Val{:diagnostic}())
    savevals = SavedValues(Float64, typeof(stateproto))
    savefunc = (u,t,integrator) -> deepcopy(getstates(setup, get_du(integrator), u, integrator.p, t, Val{:diagnostic}()))
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
	ODEProblem(func,u0,tspan,p,CryoGridODEProblem(); callback=callbacks, kwargs...)
end
# this version converts tspan from DateTime to float
"""
    CryoGridProblem(setup::CryoGridSetup, tspan::NTuple{2,DateTime}, args...;kwargs...)
"""
CryoGridProblem(setup::CryoGridSetup, tspan::NTuple{2,DateTime}, args...;kwargs...) = CryoGridProblem(setup,convert_tspan(tspan),args...;kwargs...)

export CryoGridProblem

"""
    JacobianStyle

Trait for indicating Jacobian sparsity of a CryoGrid ODEProblem.
"""
abstract type JacobianStyle end
struct DefaultJac <: JacobianStyle end
struct TridiagJac <: JacobianStyle end
"""
    JacobianStyle(::Type{<:CryoGridSetup})

Can be overriden/extended to specify Jacobian structure for specific `CryoGridSetup`s.
"""
JacobianStyle(::Type{<:CryoGridSetup}) = DefaultJac()

export DefaultJac, TridiagJac, JacobianStyle

"""
    odefunction(setup::CryoGridSetup, M, u0, p, tspan; kwargs...)

Constructs a SciML `ODEFunction` given the model setup, mass matrix M, initial state u0, parameters p, and tspan.
Can (and should) be overridden by users to provide customized ODEFunction configurations for specific problem setups, e.g:
```
model = CryoGridSetup(strat,grid)
function CryoGrid.Setup.odefunction(::DefaultJac, setup::typeof(model), M, u0, p, tspan)
    ...
    # make sure to return an instance of ODEFunction
end
...
prob = CryoGridProblem(model, tspan, p)
```

`JacobianStyle` can also be extended to create custom traits which can then be applied to compatible `CryoGridSetup`s.
"""
odefunction(setup::TSetup, M, u0, p, tspan; kwargs...) where {TSetup<:CryoGridSetup} = odefunction(JacobianStyle(TSetup), setup, M, u0, p, tspan; kwargs...)
odefunction(::DefaultJac, setup::TSetup, M, u0, p, tspan; kwargs...) where {TSetup<:CryoGridSetup} = ODEFunction(setup, mass_matrix=M; kwargs...)
function odefunction(::TridiagJac, setup::CryoGridSetup, M, u0, p, tspan; kwargs...)
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

# Auto-detect Jacobian sparsity for problems with one or more heat-only layers.
# Note: This assumes that the processes/forcings on the boundary layers do not violate the tridiagonal structure!
# Unfortunately, the Stratigraphy type signature is a bit nasty to work with :(
const HeatOnlySetup = CryoGridSetup{<:Stratigraphy{N,<:Tuple{TTop,Vararg{<:Union{<:StratNode{<:SubSurface, <:System{<:Tuple{<:Heat}}},TBot}}}}} where {N,TTop,TBot}
JacobianStyle(::Type{<:HeatOnlySetup}) = TridiagJac()
