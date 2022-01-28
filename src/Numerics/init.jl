abstract type VarInit{varname} end
Numerics.varname(::VarInit{varname}) where {varname} = varname
ConstructionBase.constructorof(::Type{T}) where {var,T<:VarInit{var}} = (args...) -> T.name.wrapper(var, args...)
"""
    ConstantInit{varname,T} <: VarInit

Initializes a scalar or vector-valued state variable with a pre-specified constant value.
"""
struct ConstantInit{varname,T} <: VarInit{varname}
    value::T
    ConstantInit(varname::Symbol, value::T) where {T} = new{varname,T}(value)
end
(init::ConstantInit)(x::AbstractVector, args...) = @. x = init.value
"""
    InterpInit{varname,P,I,E} <: VarInit{varname}

Initializer for on-grid variables that takes a `Profile` as initial values and interpolates along the model grid.
The interpolation mode is linear by default, but can also be any other `Gridded` interpolation mode supported
by `Interpolations.jl`.
"""
struct InterpInit{varname,P,I,E} <: VarInit{varname}
    profile::P
    interp::I
    extrap::E
    InterpInit(varname::Symbol, profile::P, interp::I=Linear(), extrap::E=Flat()) where {P<:Profile,I,E} = new{varname,P,I,E}(profile, interp, extrap)
end
function (init::InterpInit)(x::AbstractVector, grid::AbstractVector, args...)
    z = collect(map(knot -> dustrip(knot.depth), init.profile.knots))
    f = extrapolate(interpolate((z,), collect(map(knot -> dustrip(knot.value), init.profile.knots)), Gridded(init.interp)), init.extrap)
    @. x = f(grid)
end
# automatic partitioning of profile based on interval
Base.getindex(init::InterpInit{var}, itrv::Interval) where var = InterpInit(var, init.profile[itrv], init.interp, init.extrap)
# concatenation (violations of monotonicity are not checked...)
Base.cat(init1::InterpInit{var,P1,I,E}, init2::InterpInit{var,P2,I,E}; dims=1) where {var,P1,P2,I,E} = InterpInit(var, Profile(tuplejoin(init1.profile.knots, init2.profile.knots)), init1.interp, init1.extrap)
"""
    initializer(varname::Symbol, value::T) where {T} => ConstantInit
    initializer(varname::Symbol, profile::Profile, interp=Linear(), extrap=Flat()) => InterpInit

Convenience constructor for `VarInit`s that selects the appropriate initializer type based on the arguments.
"""
initializer(varname::Symbol, value::T) where {T} = ConstantInit(varname, value)
initializer(varname::Symbol, profile::Profile, interp=Linear(), extrap=Flat()) = InterpInit(varname, profile, interp, extrap)
