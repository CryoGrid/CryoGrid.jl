"""
    module Utils

Common utility functions, constants, and macros used throughout the CryoGrid.jl codebase.
"""
module Utils

using Dates
using ModelParameters
using StructTypes
using Unitful

import CryoGrid
import ExprTools
import ForwardDiff

export @xu_str, @Float_str, @Real_str, @Number_str, @UFloat_str, @UT_str, @setscalar, @threaded
include("macros.jl")

export DistUnit, DistQuantity, TempUnit, TempQuantity, TimeUnit, TimeQuantity
export dustrip, duconvert, applyunit
export structiterate, getscalar, tuplejoin, convert_tspan
export Params

# Convenience constants for units
const DistUnit{N} = Unitful.FreeUnits{N,Unitful.𝐋,nothing} where {N}
const DistQuantity{T,U} = Quantity{T,Unitful.𝐋,U} where {T,U<:DistUnit}
const TempUnit{N,A} = Unitful.FreeUnits{N,Unitful.𝚯,A} where {N,A}
const TempQuantity{T,U} = Quantity{T,Unitful.𝚯,U} where {T,U<:TempUnit}
const TimeUnit{N,A} = Unitful.FreeUnits{N,Unitful.𝐓,A} where {N,A}
const TimeQuantity{T,U} = Quantity{T,Unitful.𝐓,U} where {T,U<:TempUnit}

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
Base type for composite parameter types. Permits iteration of struct fields.
"""
abstract type Params end
# scalar broadcasting of Params types
Base.Broadcast.broadcastable(p::Params) = Ref(p)
# provide length and iteration over field names
Base.length(p::Params) = fieldcount(typeof(p))
Base.iterate(p::Params) = structiterate(p)
Base.iterate(p::Params, state) = structiterate(p,state)

@inline tuplejoin(x) = x
@inline tuplejoin(x, y) = (x..., y...)
@inline tuplejoin(x, y, z...) = (x..., tuplejoin(y, z...)...)

getscalar(x::Number) = x
getscalar(a::AbstractArray) = a[1]

"""
    convert_tspan(tspan::Tuple{DateTime,DateTime})
    convert_tspan(tspan::Tuple{Float64,Float64})

Convenience method for converting between `Dates.DateTime` and solver time.
"""
convert_tspan(tspan::NTuple{2,DateTime}) = Dates.datetime2epochms.(tspan) ./ 1000.0
convert_tspan(tspan::NTuple{2,Float64}) = Dates.epochms2datetime.(tspan.*1000.0)

"""
    argnames(f, choosefn=first)

Retrieves the argument names of function `f` via metaprogramming and `ExprTools`.
The optional argument `choosefn` allows for customization of which method instance
of `f` (if there is more than one) is chosen.
"""
function argnames(f, choosefn=first)
    # Parse function parameter names using ExprTools
    fms = ExprTools.methods(f)
    symbol(arg::Symbol) = arg
    symbol(expr::Expr) = expr.args[1]
    argnames = map(symbol, ExprTools.signature(choosefn(fms))[:args])
    return argnames
end

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
@generated function genmap(f, args::T) where {T<:Tuple}
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

duconvert(u::Unitful.Units, x::Number) = CryoGrid.CRYOGRID_DEBUG ? x : uconvert(u,x)

end