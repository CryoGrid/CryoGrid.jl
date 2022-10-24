"""
    Utils

Common utility functions, constants, and macros used throughout the CryoGrid.jl codebase.
"""
module Utils

using Dates
using Flatten
using IntervalSets
using ModelParameters
using Setfield
using StructTypes
using Unitful

import CryoGrid
import ForwardDiff
import FreezeCurves
import Unitful

import ConstructionBase
import FreezeCurves: normalize_temperature
import ModelParameters: stripunits

export @UFloat_str, @UT_str, @setscalar, @threaded, @sym_str, @pstrip
include("macros.jl")

export Named
export DistUnit, DistQuantity, TempUnit, TempQuantity, TimeUnit, TimeQuantity
export StrictlyPositive, StrictlyNegative, Nonnegative, Nonpositive
export dustrip, duconvert, applyunits, normalize_units, normalize_temperature, pstrip
export fastmap, fastiterate, structiterate, getscalar, tuplejoin, convert_t, convert_tspan, haskeys

# Convenience constants for units
const DistUnit{N} = Unitful.FreeUnits{N,Unitful.𝐋,nothing} where {N}
const DistQuantity{T,U} = Quantity{T,Unitful.𝐋,U} where {T,U<:DistUnit}
const TempUnit{N,A} = Unitful.FreeUnits{N,Unitful.𝚯,A} where {N,A}
const TempQuantity{T,U} = Quantity{T,Unitful.𝚯,U} where {T,U<:TempUnit}
const TimeUnit{N,A} = Unitful.FreeUnits{N,Unitful.𝐓,A} where {N,A}
const TimeQuantity{T,U} = Quantity{T,Unitful.𝐓,U} where {T,U<:TimeUnit}
# Variable/parameter domains
const StrictlyPositive = OpenInterval(0,Inf)
const StrictlyNegative = OpenInterval(-Inf,0)
const Nonnegative = Interval{:closed,:open}(0,Inf)
const Nonpositive = Interval{:open,:closed}(-Inf,0)

# StructTypes dispatches for Unitful types to convert them to and from strings
StructTypes.StructType(::Type{<:Quantity}) = StructTypes.CustomStruct()
StructTypes.lower(value::Quantity) = string(value)
StructTypes.lowertype(value::Type{<:Quantity}) = String
StructTypes.construct(::Type{Q}, value::String) where {Q<:Quantity} = uconvert(Q, uparse(replace(value, " " => "")))

"""
    Named{name,T}

Wraps an object of type `T` with a `name` type parameter.
"""
struct Named{name,T}
    obj::T
    Named(name::Symbol, obj::T) where T = new{name,T}(obj)
end
Named(values::Pair{Symbol,T}) where T = Named(values[1], values[2])
Base.nameof(::Named{name}) where name = name
ConstructionBase.constructorof(::Type{<:Named{name}}) where name = obj -> Named(name, obj)

"""
    applyunits(u::Unitful.Units, x::Number)

Conditionally applies unit `u` to `x` if and only if `x` is a unit-free quantity.
If `x` is a unitful quantity, asserts that the unit matches `u`.
"""
function applyunits(u::Unitful.Units, x::Number)
    if typeof(unit(x)) <: Unitful.FreeUnits{(), NoDims, nothing}
       return  x*u
    else
        @assert unit(x) == u "quantity $x has units incompatible with $u"
        return x
    end
end
"""
    fastmap(f::F, iter::NTuple{N,Any}...) where {F,N}

Same as `map` for `NTuple`s but with guaranteed type stability. `fastmap` is a `@generated`
function which unrolls calls to `f` into a loop-free tuple construction expression.
"""
@generated function fastmap(f::F, iters::NTuple{N,Any}...) where {F,N}
    expr = Expr(:tuple)
    for j in 1:N
        push!(expr.args, _genexpr(:(f), iters, j))
    end
    return expr
end
"""
    fastiterate(f!::F, iters::NTuple{N,Any}...) where {F,N}

Same as `fastmap` but simply invokes `f!` on each argument set without constructing a tuple.
"""
@generated function fastiterate(f!::F, iters::NTuple{N,Any}...) where {F,N}
    expr = Expr(:block)
    for j in 1:N
        push!(expr.args, _genexpr(:(f!), iters, j))
    end
    push!(expr.args, :(return nothing))
    return expr
end
_genexpr(f::Symbol, iters, j) where F = :($f($(map(i -> :(iters[$i][$j]), 1:length(iters))...)))

# special case: make sure temperatures are in °C
normalize_units(x::Unitful.AbstractQuantity{T,Unitful.𝚯}) where T = uconvert(u"°C", x)
normalize_units(x::Unitful.AbstractQuantity) = upreferred(x)
# Add method dispatch for normalize_temperature in FreezeCurves.jl
normalize_temperature(x::Param) = normalize_temperature(stripparams(x))

"""
    tuplejoin([x, y], z...)

Concatenates one or more tuples together; should generally be type stable.
"""
@inline tuplejoin(x) = x
@inline tuplejoin(x, y) = (x..., y...)
@inline tuplejoin(x, y, z...) = (x..., tuplejoin(y, z...)...)

"""
    getscalar(x)
    getscalar(x, i)

Helper method for generalizing between arrays and scalars. Without an index, retrieves
the first element of `x` if `x` is an array, otherwise simply returning `x`. If an index
`i`, is specified, returns the `i`th value of `x` if `x` is an array, or `x` otherwise.
Note that this method is not strictly necessary since Julia allows for scalar quantities
to be accessed at the first index like an array; however, the point is to make it
expliclty clear in scalar-typed code that a state variable is treated as such and is
not a vector valued quantity.
"""
getscalar(x::Number) = x
getscalar(x::Number, i) = x
getscalar(x::AbstractArray) = begin @assert length(x) == 1; x[1] end
getscalar(x::AbstractArray, i) = x[i]

"""
    convert_t(t::DateTime)
    convert_t(t::Float64)

Convenience method for converting between `Dates.DateTime` and solver time.
"""
convert_t(t::DateTime) = Dates.datetime2epochms(t) / 1000
convert_t(t::Float64) = Dates.epochms2datetime(floor(1000t))
"""
    convert_tspan(tspan::Tuple{DateTime,DateTime})
    convert_tspan(tspan::Tuple{Float64,Float64})

Convenience method for converting between `Dates.DateTime` and solver time.
"""
convert_tspan(tspan::NTuple{2,DateTime}) = convert_t.(tspan)
convert_tspan(tspan::NTuple{2,Float64}) = convert_t.(tspan)

"""
    ffill!(x::AbstractVector{T}) where {E,T<:Union{Missing,E}}

Forward fills missing values in vector `x`.
"""
function ffill!(x::AbstractVector{T}) where {E,T<:Union{Missing,E}}
    local lastval::Union{Missing,E} = missing
    @inbounds for i in 1:length(x)
        lastval = ismissing(x[i]) ? lastval : x[i]
        x[i] = lastval
    end
    return x
end

"""
Debug ustrip. Remove units if and only if debug mode is NOT enabled.
"""
dustrip(x::Number) = CryoGrid.CRYOGRID_DEBUG ? x : ustrip(x)
dustrip(x::AbstractVector{<:Quantity{T}}) where {T} = CryoGrid.CRYOGRID_DEBUG ? x : reinterpret(T, x)
dustrip(u::Unitful.Units, x::Number) = CryoGrid.CRYOGRID_DEBUG ? x : ustrip(u,x)
dustrip(u::Unitful.Units, x::AbstractVector{<:Quantity{T}}) where {T} = CryoGrid.CRYOGRID_DEBUG ? x : reinterpret(T, uconvert.(u, x))
"""
Debug uconvert.
"""
duconvert(u::Unitful.Units, x::Number) = CryoGrid.CRYOGRID_DEBUG ? x : uconvert(u, x)

"""
    pstrip(obj; keep_units=false)

Strips `Param` types and units from `obj`. If `keep_units=true`, then `Param` types will be stripped but units preserved.
"""
function pstrip(obj; keep_units=false)
    stripped_obj = ModelParameters.stripparams(obj)
    return keep_units ? stripped_obj : ModelParameters.stripunits(stripped_obj)
end

# TODO: this should be in ModelParameters.jl, not here.
function Unitful.uconvert(u::Unitful.Units, p::Param)
    nt = parent(p)
    @set! nt.val = ustrip(u, stripparams(p))
    @set! nt.units = u
    return Param(nt)
end
"""
    ModelParameters.stripunits(obj)

Additional override for `stripunits` which reconstructs `obj` with all fields that have unitful quantity
types converted to base SI units and then stripped to be unit free.
"""
function ModelParameters.stripunits(obj)
    values = Flatten.flatten(obj, Flatten.flattenable, Unitful.AbstractQuantity, Flatten.IGNORE)
    return Flatten.reconstruct(obj, map(ustrip ∘ normalize_units, values), Unitful.AbstractQuantity, Flatten.IGNORE)
end

end