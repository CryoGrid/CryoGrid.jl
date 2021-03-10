using LinearAlgebra

"""
CryoGrid specialized constructor for ODEProblem that automatically generates the initial
condition and necessary callbacks.
"""
function CryoGridProblem(setup::CryoGridSetup, tspan, p=nothing;kwargs...)
	p = isnothing(p) ? setup.pproto : p
	# compute initial condition
	u0,_ = initialcondition!(setup, p)
	N = length(u0)
	func = ODEFunction(setup;jac_prototype=Tridiagonal(ones(N-1),ones(N),ones(N-1)))
	ODEProblem(func,u0,tspan,p,kwargs...)
end

# @generated function selectvars(varnames::NamedTuple, setup::CryoGridSetup)
# 	layers = Tuple(varnames.parameters[1])
# 	select(::Type{Val{V}}) where V = (V,)
# 	select(ts::Type{<:Tuple}) = tuplejoin(map(select,Tuple(ts.parameters))...)
# 	names = (select(typeparam) for typeparam in Tuple(varnames.parameters[2].parameters))
# 	refs = (:(copy(setup.state.$layer.$name)) for (layer,names) in zip(layers,names) for name in names)
# 	:(tuple($(refs...)))
# end

export CryoGridProblem
