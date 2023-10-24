Base.@kwdef struct LinearAccumulation{S} <: SnowAccumulationScheme
    rate_scale::S = 1.0 # scaling factor for snowfall rate
end

Base.@kwdef struct DegreeDayMelt{Tfactor,Tmax} <: SnowAblationScheme
    factor::Tfactor = 5.0u"mm/K/d"
    max_unfrozen::Tmax = 0.5
end

"""
    calculate_degree_day_snow_melt(ddm::DegreeDayMelt, T_ub::Number)

Implementation of degree day melting scheme.
"""
function calculate_degree_day_snow_melt(ddm::DegreeDayMelt, T_ub::Number)
    ddf = ddm.factor # [m/K/s]
    Tref = 0.0*unit(T_ub) # just in case T_ub has units
    # calculate the melt rate per second via the degree day model
    dmelt = ddf*max(T_ub-Tref, zero(T_ub)) # [m/s]
    return max(dmelt, zero(dmelt))
end
