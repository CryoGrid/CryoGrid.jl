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
import Unitful

import ModelParameters: stripunits

export @xu_str, @Float_str, @Real_str, @Number_str, @UFloat_str, @UT_str, @setscalar, @threaded, @sym_str, @pstrip
include("macros.jl")

export DistUnit, DistQuantity, TempUnit, TempQuantity, TimeUnit, TimeQuantity
export dustrip, duconvert, applyunit, normalize_temperature, pstrip
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
    applyunit(u::Unitful.Units, x::Number)

Conditionally applies unit `u` to `x` if and only if `x` is a unit-free quantity.
If `x` is a unitful quantity, asserts that the unit matches `u`.
"""
function applyunit(u::Unitful.Units, x::Number)
    if typeof(unit(x)) <: Unitful.FreeUnits{(), NoDims, nothing}
       return  x*u
    else
        @assert unit(x) == u "quantity $x has units incompatible with $u"
        return x
    end
end

"""
    normalize_temperature(x)

Converts temperature `x` to Kelvin. If `x` has units, `uconvert` is used. Otherwise, if `x` a general numeric type, it is assumed that `x` is in celsius.
"""
normalize_temperature(x) = x + 273.15
normalize_temperature(x::TempQuantity) = uconvert(u"K", x)

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
convert_t(t::Float64) = Dates.epochms2datetime(1000t)
"""
    convert_tspan(tspan::Tuple{DateTime,DateTime})
    convert_tspan(tspan::Tuple{Float64,Float64})

Convenience method for converting between `Dates.DateTime` and solver time.
"""
convert_tspan(tspan::NTuple{2,DateTime}) = convert_t.(tspan)
convert_tspan(tspan::NTuple{2,Float64}) = convert_t.(tspan)

"""
    @generated selectat(i::Int, f, args::T) where {T<:Tuple}

Helper function for handling mixed vector/scalar arguments to runtime generated functions.
select calls getindex(i) for all array-typed arguments leaves non-array arguments as-is.
We use a generated function to expand the arguments into an explicitly defined tuple to preserve type-stability (i.e. it's an optmization);
function `f` is then applied to each element.
"""
@generated selectat(i::Int, f, args::T) where {T<:Tuple} = :(tuple($([typ <: AbstractArray ?  :(f(args[$k][i])) : :(f(args[$k])) for (k,typ) in enumerate(Tuple(T.parameters))]...)))

"""
    @generated genmap(f, args::T) where {T<:Tuple}

Generated `map` for `Tuple` types. This function is for use in generated functions where
generators/comprehensions like `map` are not allowed.
"""
@inline @generated function genmap(f, args::T) where {T<:Tuple}
    return Expr(:tuple, (:(f(args[$i])) for i in 1:length(T.parameters))...)
end

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
    adstrip(x::Number)
    adstrip(x::ForwardDiff.Dual)
    adstrip(x::ReverseDiff.TrackedReal)

adstrip extracts the underlying numeric value from `x` if `x` is a tracked value from
an autodiff library (e.g. ForwardDiff or ReverseDiff). If `x` is not an AD type, then
`adstrip` simply returns `x` as-is.
"""
adstrip(x::Number) = x
adstrip(x::ForwardDiff.Dual) = ForwardDiff.value(x) |> adstrip
adstrip(x::Param{T}) where {T<:ForwardDiff.Dual} = Param(NamedTuple{Tuple(keys(x))}((adstrip(x.val), Base.tail(parent(x))...)))

"""
Debug ustrip. Remove units if and only if debug mode is NOT enabled.
"""
dustrip(x::Number) = CryoGrid.CRYOGRID_DEBUG ? x : ustrip(x)
dustrip(x::AbstractVector{<:Quantity{T}}) where {T} = CryoGrid.CRYOGRID_DEBUG ? x : reinterpret(T, x)
dustrip(u::Unitful.Units, x::Number) = CryoGrid.CRYOGRID_DEBUG ? x : ustrip(u,x)
dustrip(u::Unitful.Units, x::AbstractVector{<:Quantity{T}}) where {T} = CryoGrid.CRYOGRID_DEBUG ? x : reinterpret(T, uconvert.(u, x))

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

Additional override for `stripunits` which reconstructs `obj` with all `AbstractQuantity` fields stripped of units.
"""
function ModelParameters.stripunits(obj)
    values = Flatten.flatten(obj, Flatten.flattenable, Unitful.AbstractQuantity, Flatten.IGNORE)
    return Flatten.reconstruct(obj, map(ustrip, values), Unitful.AbstractQuantity, Flatten.IGNORE)
end

end