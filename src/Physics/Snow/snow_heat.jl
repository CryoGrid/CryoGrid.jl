Base.@kwdef struct SturmExponential{Ta,Tb} <: SnowThermalConductivity
    a::Ta = 2.650
    b::Tb = -1.652
end

function Heat.thermalconductivity(se::SturmExponential, ρsn::Number)
    # convert ρsn to g/cm^3
    ρsn = 1000*ρsn / 1e6
    k_eff = 10^(se.a*ρsn + se.b)
    return k_eff
end

Base.@kwdef struct SturmQuadratic{Thi,Tlo,Tthresh,Tmax} <: SnowThermalConductivity
    coefs_hi::Thi = (a = 3.233, b = -1.01, c = 0.138)
    coefs_lo::Tlo = (a = 0.0, b = 0.234, c = 0.023)
    thresh::Tthresh = 0.156u"g/cm^3"
    k_max::Tmax = 2.2u"W/m/K" # assumed cond. of ice
end

function Heat.thermalconductivity(sq::SturmQuadratic, ρsn::Number)
    coefs = map((hi,lo) -> hi*(ρsn >= sq.thresh) + lo*(ρsn < sq.thresh), sq.coefs_hi, sq.coefs_lo)
    # convert ρsn to g/cm^3
    ρsn = ustrip(u"g/cm^3", applyunits(u"kg/m^3", ρsn))
    # calculate conductivity
    k_eff = coefs.c + coefs.b*ρsn + coefs.a*ρsn^2
    return min(k_eff*unit(sq.k_max), sq.k_max)
end

# extract thermal conductivity scheme from thermal properties struct and invoke special dispatches defined above.
Heat.thermalconductivity(snow::Snowpack, state, i) = thermalconductivity(snow.para.heat.cond, state.ρsn[i])

function Heat.enthalpyinv(::Snowpack, heat::HeatBalance{FreeWater,<:EnthalpyBased}, H, θwi, C, L)
    T_f = H / C
    T = IfElse.ifelse(
        H < zero(θwi),
        # Case 1: H < 0 -> frozen
        T_f,
        # Case 2: H >= 0
        # For snow, the temperature is not permitted to rise above zero
        # since H > 0 implies snow melt. In cases where degree day melt
        # is used, this may lead to some excess energy in the snowpack.
        zero(T_f),
    )
    return T
end

Heat.freezethaw!(snowpack::Snowpack, state) = Heat.freezethaw!(snowpack, snowpack.heat, state)
