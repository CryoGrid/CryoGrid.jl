abstract type VarInitializer{varname} end
Numerics.varname(::VarInitializer{varname}) where {varname} = varname
ConstructionBase.constructorof(::Type{T}) where {varname,T<:VarInitializer{varname}} = (args...) -> T.name.wrapper(varname, args...)
"""
    ConstantInitializer{varname,T} <: VarInitializer{varname}

Initializes a scalar or vector-valued state variable with a pre-specified constant value.
"""
struct ConstantInitializer{varname,T} <: VarInitializer{varname}
    value::T
    ConstantInitializer(varname::Symbol, value::T) where {T} = new{varname,T}(value)
end
(init::ConstantInitializer)(u::AbstractVector, args...) = @. u = init.value
Base.getindex(init::ConstantInitializer, itrv::Interval) = init
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
function (init::InterpInitializer)(u::AbstractVector, z::AbstractVector)
    @unpack profile, interp, extrap = init
    depths = collect(map(knot -> dustrip(knot.depth), profile.knots))
    if length(depths) > 1
        f = extrapolate(interpolate((depths,), collect(map(knot -> dustrip(knot.value), profile.knots)), Gridded(interp)), extrap)
        @. u = f(z)
        return u
    else
        # if only one knot is defined, set to this value over all z;
        # this is necessary because interpolate does not support interpolation grids of length 1
        y = dustrip(profile.knots[1].value)
        @. u = y
        return u
    end
end
# automatic partitioning of profile based on interval
Base.getindex(init::InterpInitializer{var}, itrv::Interval) where var = InterpInitializer(var, init.profile[itrv], init.interp, init.extrap)
"""
    initializer(varname::Symbol, value::T) where {T} => ConstantInitializer
    initializer(varname::Symbol, profile::Profile, interp=Linear(), extrap=Flat()) => InterpInitializer

Convenience constructor for `VarInitializer` that selects the appropriate initializer type based on the arguments.
"""
initializer(varname::Symbol, value::T) where {T} = ConstantInitializer(varname, value)
initializer(varname::Symbol, profile::Profile, interp=Linear(), extrap=Flat()) = InterpInitializer(varname, profile, interp, extrap)
