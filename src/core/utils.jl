# Convenience constants for units
const Unit{N,D,A} = Unitful.Units{N,D,A}
const DistUnit{N} = Unitful.FreeUnits{N,Unitful.𝐋,nothing} where {N}
const DistQuantity{T,U} = Quantity{T,Unitful.𝐋,U} where {T,U<:DistUnit}
const TempUnit{N,A} = Unitful.FreeUnits{N,Unitful.𝚯,A} where {N,A}
const TempQuantity{T,U} = Quantity{T,Unitful.𝚯,U} where {T,U<:TempUnit}
const TimeUnit{N,A} = Unitful.FreeUnits{N,Unitful.𝐓,A} where {N,A}
const TimeQuantity{T,U} = Quantity{T,Unitful.𝐓,U} where {T,U<:TempUnit}

export Unit, DistUnit, DistQuantity, TempUnit, TempQuantity

"""
Similar to Unitful.@u_str (i.e. u"kg") but conditional on debug mode being enabled. Otherwise, no unit is applied.
This should be used to apply units (and thus dimensional analysis checks) to physical quantities at test time but
not during normal execution to avoid unnecessary overhead.
"""
macro xu_str(unit) CRYOGRID_DEBUG ? :(@u_str($unit)) : 1 end
"""
Similar to @UT_str but produces a Float64 quantity type for the given unit if and only if debug mode is enabled.
If debug mode is not enabled, plain Float64 is used instead.
"""
macro Float_str(unit) CRYOGRID_DEBUG ? :(typeof(@u_str($unit)*0.0)) : :(Float64) end
"""
Similar to @UT_str but produces a Real quantity type for the given unit if and only if debug mode is enabled.
If debug mode is not enabled, plain Float64 is used instead.
"""
macro Real_str(unit) CRYOGRID_DEBUG ? :(typeof(@u_str($unit)*0.0)) : :(Real) end
"""
Similar to Unitful.@u_str (i.e. u"kg") but produces the type of the quantity rather than the instance. NOT conditional
on debug mode.
"""
macro UFloat_str(unit) :(typeof(@u_str($unit)*0.0)) end
"""
Similar to Unitful.@u_str (i.e. u"kg") but produces the type of the unit rather than the instance. NOT conditional
on debug mode.
"""
macro UT_str(unit) :(typeof(@u_str($unit))) end

"""
Debug ustrip. Remove units if and only if debug mode is NOT enabled.
"""
dustrip(x::Number) = CRYOGRID_DEBUG ? x : ustrip(x)
dustrip(u::Unitful.Units, x::Number) = CRYOGRID_DEBUG ? x : ustrip(u,x)

export @xu_str, @Float_str, @UFloat_str, @UT_str, dustrip

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

export structiterate

@inline tuplejoin(x) = x
@inline tuplejoin(x, y) = (x..., y...)
@inline tuplejoin(x, y, z...) = (x..., tuplejoin(y, z...)...)

export tuplejoin

"""
Convenience macro for setting scalar (single-element) arrays/vectors. It turns an expression of the form:
    `a.b = ...`
into
    `a.b[1] = ...`

This is primarily intended for code clarity, i.e to clearly discern scalar and non-scalar values.
"""
macro setscalar(expr)
    refexpr = expr.args[1]
    valexpr = expr.args[2]
    quote
        $(esc(refexpr))[1] = $(esc(valexpr))
    end
end

getscalar(x::Number) = x
getscalar(a::AbstractArray) = a[1]

export @setscalar, getscalar

"""
    Profile(pairs...;names)

Constructs a Profile from the given pairs Q => (x1,...,xn) where x1...xn are the values defined at Q.
Column names for the resulting DimArray can be set via the names parameter which accepts an NTuple of symbols,
where N must match the number of parameters given (i.e. n).
"""
function Profile(pairs::Pair{Q,NTuple{N,T}}...;names::Union{Nothing,NTuple{N,Symbol}}=nothing) where {T,N,Q<:DistQuantity}
    depths, vals = zip(pairs...)
    params = hcat(collect.(vals)...)'
    names = isnothing(names) ? [Symbol(:x,:($i)) for i in 1:N] : collect(names)
    DimArray(params, (Z(collect(depths)), Y(names)))
end

"""
    interpolateprofile(profile::Profile, state; interp=Linear())

Interpolates the given profile to the corresponding variable grids. Assumes state to be indexable via the corresponding
variable symbol and that the parameter names in state and profile match.
"""
function interpolateprofile!(profile::DimArray, state; interp=Linear())
    let (depths,names) = dims(profile),
        z = ustrip.(depths);
        for p in names
            # in case state is unit-free, reinterpret to match eltype of profile
            pstate = reinterpret(eltype(profile),state[p])
            f = @> interpolate((z,), profile[:,p], Gridded(interp)) extrapolate(Flat())
            pstate .= f.(state.grids[p])   # assume length(grid) == length(state.p)
        end
    end
end

export Profile, interpolateprofile!

"""
    convert_tspan(tspan::Tuple{DateTime,DateTime})
    convert_tspan(tspan::Tuple{Float64,Float64})

Convenience method for converting between `Dates.DateTime` and solver time.
"""
convert_tspan(tspan::NTuple{2,DateTime}) = Dates.datetime2epochms.(tspan) ./ 1000.0
convert_tspan(tspan::NTuple{2,Float64}) = Dates.epochms2datetime.(tspan.*1000.0)
export convert_tspan

# Helper function for handling mixed vector/scalar arguments to runtime generated functions.
# select calls getindex(i) for all array-typed arguments leaves non-array arguments as-is.
# We use a generated function to expand the arguments into an explicitly defined tuple to preserve type-stability (i.e. it's an optmization);
# function f is then applied to each element
@generated selectat(i::Int, f, args::T) where {T<:Tuple} = :(tuple($([typ <: AbstractArray ?  :(f(args[$k][i])) : :(f(args[$k])) for (k,typ) in enumerate(Tuple(T.parameters))]...)))

"""
adstrip extracts the underlying numeric value from `x` if `x` is a tracked value from
an autodiff library (e.g. ForwardDiff or ReverseDiff). If `x` is not an AD type, then
`adstrip` simply returns `x` as-is.
"""
adstrip(x::Number) = x
adstrip(x::ForwardDiff.Dual) = ForwardDiff.value(x)
adstrip(x::ReverseDiff.TrackedReal) = x.value
