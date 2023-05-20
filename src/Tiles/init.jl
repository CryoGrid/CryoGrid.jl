# default behavior is to not automatically parameterize initializers
CryoGrid.parameterize(init::VarInitializer) = init

"""
    FunctionInitializer{varname,F} <: VarInitializer{varname}

Initializes a scalar or vector-valued state variable using an arbitrary method/function.
"""
struct FunctionInitializer{varname,F} <: VarInitializer{varname}
    f::F
    FunctionInitializer(varname::Symbol, f::F) where {F} = new{varname,F}(f)
end
CryoGrid.initialcondition!(layer::Layer, process::Process, state, init::FunctionInitializer) = init.f(layer, process, state)
Base.getindex(init::FunctionInitializer, itrv::Interval) = init

"""
    InterpInitializer{varname,P,I,E} <: VarInitializer{varname}

Initializer for on-grid variables that takes a `Profile` as initial values and interpolates along the model grid.
The interpolation mode is linear by default, but can also be any other `Gridded` interpolation mode supported
by `Interpolations.jl`.
"""
struct InterpInitializer{varname,P,I,E} <: VarInitializer{varname}
    profile::P
    interp::I
    extrap::E
    InterpInitializer(varname::Symbol, profile::P, interp::I=Linear(), extrap::E=Flat()) where {P<:Profile,I,E} = new{varname,P,I,E}(profile, interp, extrap)
end

function CryoGrid.initialcondition!(layer::Layer, process::Process, state, init::InterpInitializer{var}) where var
    profile, interp, extrap = init.profile, init.interp, init.extrap
    depths = collect(map(knot -> ustrip(knot.depth), profile.knots))
    u = getproperty(state, var)
    z = cells(state.grid)
    @assert length(z) == length(u)
    if length(depths) > 1
        f = Interpolations.extrapolate(
            Interpolations.interpolate(
                (depths,),
                collect(map(knot -> ustrip(knot.value), profile.knots)),
                Interpolations.Gridded(interp)
            ),
            extrap
        )
        @. u = f(z)
        return u
    else
        # if only one knot is defined, set to this value over all z;
        # this is necessary because interpolate does not support interpolation grids of length 1
        y = ustrip(profile.knots[1].value)
        @. u = y
        return u
    end
end

# automatic partitioning of profile based on interval
Base.getindex(init::InterpInitializer{var}, itrv::Interval) where var = InterpInitializer(var, init.profile[itrv], init.interp, init.extrap)

"""
    initializer(varname::Symbol, x::Number) => FunctionInitializer w/ constant
    initializer(varname::Symbol, f::Function) => FunctionInitializer
    initializer(varname::Symbol, profile::Profile, interp=Linear(), extrap=Flat()) => InterpInitializer

Convenience constructor for `VarInitializer` that selects the appropriate initializer type based on the arguments.
"""
initializer(varname::Symbol, x::Number) = FunctionInitializer(varname, (layer,proc,state) -> getproperty(state, varname) .= x)
initializer(varname::Symbol, f::Function) = FunctionInitializer(varname, f)
initializer(varname::Symbol, profile::Profile, interp=Interpolations.Linear(), extrap=Interpolations.Flat()) = InterpInitializer(varname, profile, interp, extrap)
