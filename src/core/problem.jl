"""
CryoGrid specialized constructor for ODEProblem that automatically generates the initial
condition and necessary callbacks.
"""
function CryoGridProblem(setup::CryoGridSetup, tspan;savevars::NamedTuple=NamedTuple(),kwargs...)
	# "lifts" a symbol or tuple of symbols into Val types to be interpreted by the compiler
	lift(s::Symbol) = Val{s}()
	lift(ss::Tuple{Vararg{Symbol}}) = map(lift, ss)
	varvals = map(lift, savevars)
	vars = selectvars(varvals,setup)
	# set up saving callback with given variables
	cache = SavedValues(eltype(tspan), typeof(vars))
	savecallback = SavingCallback((u,t,integrator)->selectvars(varvals,setup), cache;kwargs...)
	callbacks = CallbackSet(savecallback)
	# compute initial condition
	u0,_ = initialcondition!(setup)
	ODEProblem(setup,u0,tspan,callback=callbacks), cache
end

@generated function selectvars(varnames::NamedTuple, setup::CryoGridSetup)
	layers = Tuple(varnames.parameters[1])
	select(::Type{Val{V}}) where V = (V,)
	select(ts::Type{<:Tuple}) = tuplejoin(map(select,Tuple(ts.parameters))...)
	names = (select(typeparam) for typeparam in Tuple(varnames.parameters[2].parameters))
	refs = (:(copy(setup.state.$layer.$name)) for (layer,names) in zip(layers,names) for name in names)
	:(tuple($(refs...)))
end

export CryoGridProblem
