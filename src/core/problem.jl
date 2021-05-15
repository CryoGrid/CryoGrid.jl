"""
Specialized problem type for CryoGrid ODEProblems.
"""
struct CryoGridProblem end

"""
    CryoGridProblem(setup::CryoGridSetup, tspan::NTuple{2,Float64}, p=nothing;kwargs...)

CryoGrid specialized constructor for ODEProblem that automatically generates the initial
condition and necessary callbacks.
"""
function CryoGridProblem(setup::CryoGridSetup, tspan::NTuple{2,Float64}, p=nothing;kwargs...)
	p = isnothing(p) ? setup.pproto : p
	# compute initial condition
	u0,_ = initialcondition!(setup, p, tspan)
    # set up default mass matrix
    M_diag = similar(setup.uproto)
    for layer in keys(setup.meta)
        progvars = setup.meta[layer][:progvars]
        algvars = setup.meta[layer][:algvars]
        for var in progvars
            M_diag[layer][varname(var)] .= one(eltype(M_diag))
        end
        for var in algvars
            M_diag[layer][varname(var)] .= zero(eltype(M_diag))
        end
    end
    # if no algebraic variables are present, use identity matrix
    num_algebraic = length(M_diag) - sum(M_diag)
    M = num_algebraic > 0 ? Diagonal(M_diag) : I
	func = odefunction(setup, M, u0, p, tspan; kwargs...)
	ODEProblem(func,u0,tspan,p,CryoGridProblem(); kwargs...)
end
# this version converts tspan from DateTime to float
"""
    CryoGridProblem(setup::CryoGridSetup, tspan::NTuple{2,DateTime}, args...;kwargs...)
"""
CryoGridProblem(setup::CryoGridSetup, tspan::NTuple{2,DateTime}, args...;kwargs...) = CryoGridProblem(setup,convert_tspan(tspan),args...;kwargs...)

export CryoGridProblem

"""
Trait for determining Jacobian sparsity of a CryoGrid ODEProblem.
"""
abstract type JacobianStyle end
struct DefaultJac <: JacobianStyle end
struct TridiagJac <: JacobianStyle end
# struct SparseJac <: JacobianStyle end
JacobianStyle(::Type{<:CryoGridSetup}) = DefaultJac()

"""
    odefunction(setup::CryoGridSetup, M, u0, p, tspan; kwargs...)

Constructs a SciML `ODEFunction` given the model setup, mass matrix M, initial state u0, parameters p, and tspan.
Can (and should) be overridden by users to provide customized ODEFunction configurations for specific problem setups, e.g:
```
model = CryoGridSetup(strat,grid)
function CryoGrid.odefunction(::DefaultJac, setup::typeof(model), M, u0, p, tspan)
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
        ODEFunction(setup, mass_matrix=M, kwargs...)
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
