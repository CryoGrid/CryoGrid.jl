using CryoGrid
using CryoGrid: Parameterization, DynamicParameterization
using Interpolations
using ModelParameters

import ConstructionBase

# Trend

"""
    LinearTrend{id,TSlope,TInt} <: DynamicParameterization

Applies a linear trend to a parameter `p` by reparameterizing it as: `p = p₁*t + p₀`
"""
struct LinearTrend{id,TSlope,TInt} <: DynamicParameterization
    slope::TSlope
    intercept::TInt
    tstart::Float64
    tstop::Float64 
    period::Float64
    minval::Float64
    maxval::Float64
    function LinearTrend(id::Symbol, slope=Param(0.0), intercept=Param(0.0); tstart=0.0, tstop=Inf, period=1.0, minval=-Inf, maxval=Inf)
        return new{id,typeof(slope),typeof(intercept)}(slope, intercept, tstart, tstop, period, minval, maxval)
    end
end

function (trend::LinearTrend)(t)
    let t = min(max(t - trend.tstart, zero(t)), trend.tstop - trend.tstart),
        β = trend.slope / trend.period,
        α = trend.intercept;
        min(max(β*t + α, trend.minval), trend.maxval)
    end
end


ModelParameters.component(::Type{<:LinearTrend{id}}) where {id} = LinearTrend{id}

CryoGrid.InputOutput.paramname(::Param, ::Type{<:LinearTrend{id}}, fieldname::Symbol) where id = Symbol(id, :(.), fieldname)

# constructor overload to allow automatic reconstruction
ConstructionBase.constructorof(::Type{<:LinearTrend{id}}) where id = (slope, intercept, tstart, tstop, period, minval, maxval) -> LinearTrend(id, slope, intercept; tstart, tstop, period, minval, maxval)

# Piecewise linear

"""
Helper type for PiecewiseLinear.
"""
struct PiecewiseKnot{id,Tw,Tv}
    binwidth::Tw
    value::Tv
    PiecewiseKnot(id::Symbol, binwidth, value) = new{id,typeof(binwidth),typeof(value)}(binwidth, value)
end
ConstructionBase.constructorof(::Type{<:PiecewiseKnot{id}}) where id = (w,v) -> PiecewiseKnot(id, w, v)

ModelParameters.component(::Type{<:PiecewiseKnot{id}}) where {id} = PiecewiseKnot{id}

CryoGrid.InputOutput.paramname(::Param, ::Type{<:PiecewiseKnot{id}}, fieldname::Symbol) where id = Symbol(id, :(.), fieldname)

"""
    PiecewiseLinear{id,N,T0,TKnots,I} <: DynamicParameterization

Reparameterizes parameter `p` as `p = p₁δ₁t + ⋯ + pₖδₖt` where δₖ are indicators
for when `tₖ₋₁ <= t <= tₖ`. To facilitate sampling and optimization, change points
tᵢ are parameterized as bin widths, which should be strictly positive. `PiecewiseLinear`
will normalize them and scale by the size of the time interval.
"""
struct PiecewiseLinear{id,N,T0,TKnots,I<:Union{Constant,Linear}} <: DynamicParameterization
    initialvalue::T0
    knots::TKnots
    tstart::Float64
    tstop::Float64
    interp::I
    function PiecewiseLinear(id::Symbol, initialvalue::T0, knots::TKnots, tstart::Float64, tstop::Float64, interp::I=Linear()) where {T0,TKnots,I}
        @assert tstop >= tstart
        @assert isfinite(tstart) && isfinite(tstop)
        new{id,length(knots),T0,TKnots,I}(initialvalue, knots, tstart, tstop, interp)
    end
end
PiecewiseLinear(id::Symbol, initialvalue; tstart=0.0, tstop=1.0, interp=Linear()) = PiecewiseLinear(id, initialvalue, (), tstart, tstop, interp)
PiecewiseLinear(id::Symbol, initialvalue, knots::NTuple{2,Any}...; tstart=0.0, tstop=1.0, interp=Linear()) = PiecewiseLinear(id, initialvalue, Tuple(map(knot -> PiecewiseKnot(id, knot...), knots)), tstart, tstop, interp)
"""
Special case of `PiecewiseLinear`.
"""
PiecewiseConstant(id::Symbol, args...; kwargs...) = PiecewiseLinear(id, args...; interp=Constant(), kwargs...)

ModelParameters.component(::Type{<:PiecewiseLinear{id}}) where {id} = PiecewiseLinear{id}

CryoGrid.InputOutput.paramname(::Param, ::Type{<:PiecewiseLinear{id}}, fieldname::Symbol) where id = Symbol(id, :(.), fieldname)

ConstructionBase.constructorof(::Type{<:PiecewiseLinear{id}}) where id = (i,k,t0,t1,interp) -> PiecewiseLinear(id, i, k, t0, t1, interp)

function (pc::PiecewiseLinear{id,N,T0,TKnots,I})(t) where {id,N,T0,TKnots,I<:Union{Constant,Linear}}
    f(::Type{<:Constant}, t, tᵢ, xᵢ, tᵢ₊₁, xᵢ₊₁) = xᵢ₊₁
    function f(::Type{<:Linear}, t, tᵢ, xᵢ, tᵢ₊₁, xᵢ₊₁)
        if t <= tᵢ
            return xᵢ
        elseif t >= tᵢ₊₁
            return xᵢ₊₁
        else
            return xᵢ + (xᵢ₊₁ - xᵢ)*(t - tᵢ) / (tᵢ₊₁ - tᵢ)
        end
    end
    let tspan = pc.tstop - pc.tstart,
        t = t - pc.tstart;
        if t < zero(t)
            return pc.initialvalue
        else
            bins = map(k -> k.binwidth, pc.knots)
            vals = (pc.initialvalue, map(k -> k.value, pc.knots)...)
            ts = (0.0, cumsum((bins ./ sum(bins)).*tspan)...)
            i = binindex(ts, 1, length(ts), t)
            return f(I, t, ts[i], vals[i], ts[i+1], vals[i+1])
        end
    end
end
function binindex(values::Tuple, st, en, x)
    mid = Int(floor((st + en)/2))
    if values[mid] >= x && en - st <= 1
        return st
    elseif values[mid] < x && en - st <= 1
        return en-1
    elseif values[mid] > x
        return binindex(values, st, mid, x)
    else
        return binindex(values, mid, en, x)
    end
end

# Transformed parameterization

"""
    Transformed{F,NT} <: DynamicParameterization

Reparameterization corresponding to an arbitrary functional transformation.
"""
struct Transformed{F,NT} <: DynamicParameterization
    f::F
    args::NT
    Transformed(f::F, args::NT) where {F,NT<:NamedTuple} = new{F,NT}(f, args)
    function Transformed(f::F, p::Pair{Symbol,<:Param}) where {F}
        nt = NamedTuple(tuple(p))
        new{F,typeof(nt)}(f, nt)
    end
end

(tf::Transformed)(t) = tf.f(map(stripparams, tf.args)...)

# "Time-varying" parameterization

"""
    TimeVarying{F,NT} <: DynamicParameterization

Reparameterization corresponding to an arbitrary time-varying function
with optional keyword arguments.
"""
struct TimeVarying{F,NT} <: DynamicParameterization
    f::F
    kwargs::NT
    TimeVarying(f::F, kwargs::NT=(;)) where {F,NT<:NamedTuple} = new{F,NT}(f, kwargs)
end

(tf::TimeVarying)(t) = tf.f(t; tf.kwargs...)

