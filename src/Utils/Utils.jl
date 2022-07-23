"""
    Utils

Common utility functions, constants, and macros used throughout the CryoGrid.jl codebase.
"""
module Utils

using Dates
using Flatten
using ModelParameters
using Setfield
using StructTypes
using Unitful

import CryoGrid
import ForwardDiff
import FreezeCurves
import Unitful

import ModelParameters: stripunits

export @xu_str, @Float_str, @Real_str, @Number_str, @UFloat_str, @UT_str, @setscalar, @threaded, @sym_str, @pstrip
include("macros.jl")

export DistUnit, DistQuantity, TempUnit, TempQuantity, TimeUnit, TimeQuantity
export dustrip, duconvert, applyunits, normalize_temperature, pstrip
export structiterate, getscalar, tuplejoin, convert_t, convert_tspan, haskeys
export IterableStruct

# Convenience constants for units
const DistUnit{N} = Unitful.FreeUnits{N,Unitful.𝐋,nothing} where {N}
const DistQuantity{T,U} = Quantity{T,Unitful.𝐋,U} where {T,U<:DistUnit}
const TempUnit{N,A} = Unitful.FreeUnits{N,Unitful.𝚯,A} where {N,A}
const TempQuantity{T,U} = Quantity{T,Unitful.𝚯,U} where {T,U<:TempUnit}
const TimeUnit{N,A} = Unitful.FreeUnits{N,Unitful.𝐓,A} where {N,A}
const TimeQuantity{T,U} = Quantity{T,Unitful.𝐓,U} where {T,U<:TimeUnit}

StructTypes.StructType(::Type{<:Quantity}) = StructTypes.CustomStruct()
StructTypes.lower(value::Quantity) = string(value)
StructTypes.lowertype(value::Type{<:Quantity}) = String
StructTypes.construct(::Type{Q}, value::String) where {Q<:Quantity} = uconvert(Q, uparse(replace(value, " " => "")))

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

FreezeCurves.normalize_temperature(x::Param) = stripparams(x) |> FreezeCurve.normalize_temperature

"""
Provides implementation of `Base.iterate` for structs.
"""
function structiterate(obj::A) where {A}
    names = fieldnames(A)
    if length(names) == 0; return nothing end
    gen = (getfield(obj,name) for name in names)
    (val,state) = iterate(gen)
    (val, (gen,state))
end
function structiterate(obj, state)
    gen, genstate = state
    nextitr = iterate(gen,genstate)
    isnothing(nextitr) ? nothing : (nextitr[1],(gen,nextitr[2]))
end

"""
Base type for allowing iteration of struct fields.
"""
abstract type IterableStruct end
# scalar broadcasting of IterableStruct types
Base.Broadcast.broadcastable(p::IterableStruct) = Ref(p)
# provide length and iteration over field names
Base.length(p::IterableStruct) = fieldcount(typeof(p))
Base.iterate(p::IterableStruct) = structiterate(p)
Base.iterate(p::IterableStruct, state) = structiterate(p,state)

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

duconvert(u::Unitful.Units, x::Number) = CryoGrid.CRYOGRID_DEBUG ? x : uconvert(u, x)

param(x::Number, kwargs...) = Param(x, kwargs...)
param(x::Unitful.AbstractQuantity, kwargs...) = Param(ustrip(x), units=units(x), kwargs...)
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
    # special case: make sure temperatures are in °C
    normalize_units(x::Unitful.AbstractQuantity{T,Unitful.𝚯}) where T = uconvert(u"°C", x)
    normalize_units(x::Unitful.AbstractQuantity) = upreferred(x)
    values = Flatten.flatten(obj, Flatten.flattenable, Unitful.AbstractQuantity, Flatten.IGNORE)
    return Flatten.reconstruct(obj, map(ustrip ∘ normalize_units, values), Unitful.AbstractQuantity, Flatten.IGNORE)
end

end