# Generic step limiter types
abstract type StepLimiter end
CryoGrid.parameterize(limiter::StepLimiter) = limiter

"""
    MaxDelta{T}

Allow a maximum change of `Δmax` in the integrated quantity.
"""
Base.@kwdef struct MaxDelta{T} <: StepLimiter
    Δmax::T
    upper_limit_factor::Float64 = 0.99
    lower_limit_factor::Float64 = 0.99
end
MaxDelta(Δmax) = MaxDelta(;Δmax)
function (limiter::MaxDelta)(du, u, t)
    dtmax = abs(limiter.Δmax / du)
    return isfinite(dtmax) ? dtmax : Inf
end
function (limiter::MaxDelta)(
    du, u, t, lower_limit, upper_limit;
    lower_limit_factor=limiter.lower_limit_factor,
    upper_limit_factor=limiter.upper_limit_factor,
)
    dtmax = IfElse.ifelse(
        sign(du) > 0,
        upper_limit_factor*(upper_limit - u) / du,
        -lower_limit_factor*(u - lower_limit) / du,
    )
    dtmax = min(dtmax, abs(limiter.Δmax / du))
    return isfinite(dtmax) ? dtmax : Inf
end
"""
    CFL{Tmax<:MaxDelta}

Courant-Fredrichs-Lewy condition (where defined) with the given Courant number and embedded `MaxDelta` condition.
"""
Base.@kwdef struct CFL{Tmax<:MaxDelta} <: StepLimiter
    courant_number::Float64 = 0.5
    maxdelta::Tmax = MaxDelta(Inf)
end
