abstract type VarInit{varname} end
Numerics.varname(::VarInit{varname}) where {varname} = varname
"""
    init!(x::AbstractVector, init::VarInit, args...)

Initializes state variable `x` using initializer `init` and implementation-specific additional
arguments.
"""
init!(x::AbstractVector, init::VarInit, args...) = error("not implemented")
"""
    ConstantInit{varname,T} <: VarInit

Initializes a scalar or vector-valued state variable with a pre-specified constant value.
"""
struct ConstantInit{varname,T} <: VarInit{varname}
    value::T
    ConstantInit(varname::Symbol, value::T) where {T} = new{varname,T}(value)
end
init!(x::AbstractVector, init::ConstantInit) = @. x = init.value
"""
    InterpInit{varname,P,I,E} <: VarInit

Init for on-grid variables that takes a `Profile` as initial values
and interpolates along the model grid. The interpolation mode is linear by default,
but can also be quadratic, cubic, or any other `Gridded` interpolation mode supported
by `Interpolations.jl`.
"""
struct InterpInit{varname,P,I,E} <: VarInit{varname}
    profile::P
    interp::I
    extrap::E
    InterpInit(varname::Symbol, profile::P, interp::I=Linear(), extrap::E=Flat()) where {P<:Profile,I,E} = new{varname,P,I,E}(profile, interp)
end
function init!(x::AbstractVector, init::InterpInit, grid::AbstractVector)
    z = collect(map(dustrip, init.profile.depths))
    f = extrapolate(interpolate((z,), collect(map(dustrip, init.profile.values)), Gridded(init.interp)), init.extrap)
    @. x = f(grid)
end
"""
    initializer(varname::Symbol, value::T) where {T} => ConstantInit
    initializer(varname::Symbol, profile::Profile, interp=Linear(), extrap=Flat()) => InterpInit

Convenience constructor for `VarInit`s that selects the appropriate initializer type based on the arguments.
"""
initializer(varname::Symbol, value::T) where {T} = ConstantInit(varname, value)
initializer(varname::Symbol, profile::Profile, interp=Linear(), extrap=Flat()) = InterpInit(varname, profile, interp, extrap)
