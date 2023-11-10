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

Base.@kwdef struct SturmQuadratic{Trh,Thi,Tlo} <: SnowThermalConductivity
    range_hi::Trh = 156..600.0
    coefs_hi::Thi = (a = 3.233, b = -1.01, c = 0.138)
    coefs_lo::Tlo = (a = 0.0, b = 0.234, c = 0.023)
end

function Heat.thermalconductivity(sq::SturmQuadratic, ρsn::Number)
    coefs = if ρsn ∈ sq.range_hi
        sq.coefs_hi
    elseif ρsn < infimum(sq.range_hi)
        sq.coefs_lo
    else
        error("snow density $ρsn outside of supported range (max: $(supremum(sq.range_hi)))")
    end
    # convert ρsn to g/cm^3
    ρsn = 1000*ρsn / 1e6
    # calculate conductivity
    k_eff = coefs.c + coefs.b*ρsn + coefs.a*ρsn^2
    return k_eff
end

# extract thermal conductivity scheme from thermal properties struct and invoke special dispatches defined above.
Heat.thermalconductivity(snow::Snowpack, heat::HeatBalance, state, i) = thermalconductivity(snow.para.heat.cond, state.ρsn[i])

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
