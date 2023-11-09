Base.@kwdef struct SturmExponential{Ta,Tb} <: SnowThermalConductivity
    a::Ta = 2.650
    b::Tb = -1.652
end

Heat.thermalconductivity(se::SturmExponential, ρsn::Number) = 10^(se.a*ρsn + se.b)

Base.@kwdef struct SturmQuadratic{Trh,Thi,Tlo} <: SnowThermalConductivity
    range_hi::Trh = 0.156..0.6
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
    # calculate conductivity
    k_eff = coefs.c + coefs.b*ρsn + coefs.a*ρsn^2
    return k_eff
end

# extract thermal conductivity scheme from thermal properties struct and invoke special dispatches defined above.
Heat.thermalconductivity(snow::Snowpack, heat::HeatBalance, state, i) = thermalconductivity(snow.para.heat, state.ρsn[i])
