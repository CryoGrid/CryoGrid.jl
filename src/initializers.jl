Base.isless(init1::Initializer, init2::Initializer) = false

# add varname dispatch for initializer types
CryoGrid.varname(::VarInitializer{varname}) where {varname} = varname

# default behavior is to not automatically parameterize initializers
CryoGrid.parameterize(init::VarInitializer) = init

# default to invoking initializer
CryoGrid.initialcondition!(init::VarInitializer, layer::Layer, state) = init(layer, state)

"""
    FunctionInitializer{varname,F} <: VarInitializer{varname}

Initializes a scalar or vector-valued state variable using an arbitrary method/function.
"""
struct FunctionInitializer{varname,F} <: VarInitializer{varname}
    f::F
    FunctionInitializer(varname::Symbol, f::F) where {F} = new{varname,F}(f)
end

(init::FunctionInitializer)(layer::Layer, state) = init.f(layer, state)

# do not deconstruct FunctionInitializer
Flatten.flattenable(::Type{<:FunctionInitializer}, ::Type{Val{:f}}) = false

ConstructionBase.constructorof(::Type{T}) where {varname,T<:FunctionInitializer{varname}} = f -> FunctionInitializer(varname, f)

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

function (init::InterpInitializer{var})(::Layer, state) where {var}
    profile, interp, extrap = init.profile, init.interp, init.extrap
    depths = collect(ustrip.(keys(profile)))
    u = getproperty(state, var)
    z = cells(state.grid)
    @assert length(z) == length(u) "Length of grid does not match state of $var on layer $(state.name): $(length(z)) != $(length(u))"
    if length(depths) > 1
        f = Interpolations.extrapolate(
            Interpolations.interpolate(
                (depths,),
                collect(ustrip.(values(profile))),
                Interpolations.Gridded(interp)
            ),
            extrap
        )
        f_adapt = adapt(Numerics.arraytype(state.grid), f)
        @. u = f_adapt(z)
        return u
    else
        # if only one knot is defined, set to this value over all z;
        # this is necessary because interpolate does not support interpolation grids of length 1
        y = ustrip(first(values(profile)))
        @. u = y
        return u
    end
end

ConstructionBase.constructorof(::Type{T}) where {varname,T<:InterpInitializer{varname}} = (profile, interp, extrap) -> InterpInitializer(varname, profile, interp, extrap)

# automatic partitioning of profile based on interval
Base.getindex(init::InterpInitializer{var}, itrv::Interval) where var = InterpInitializer(var, init.profile[itrv], init.interp, init.extrap)

"""
    initializer(varname::Symbol, args...) = initializer(Val{varname}(), args...)
    initializer(::Val{varname}, x::Number) => FunctionInitializer w/ constant
    initializer(::Val{varname}, f::Function) => FunctionInitializer
    initializer(::Val{varname}, profile::Profile, interp=Linear(), extrap=Flat()) => InterpInitializer

Convenience constructor for `VarInitializer` that selects the appropriate initializer type based on the arguments.
"""
initializer(varname::Symbol, args...) = initializer(Val{varname}(), args...)
initializer(::Val{varname}, x::Number) where {varname} = FunctionInitializer(varname, (layer,state) -> getproperty(state, varname) .= x)
initializer(::Val{varname}, f::Function) where {varname} = FunctionInitializer(varname, f)
initializer(::Val{varname}, profile::Profile, interp=Interpolations.Linear(), extrap=Interpolations.Flat()) where {varname} = InterpInitializer(varname, profile, interp, extrap)
# if `initializer` is passed a matching initializer, just return it as-is.
initializer(::Val{varname}, init::VarInitializer{varname}) where {varname} = init
