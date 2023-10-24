abstract type SnowDensityScheme end

# constant density (using Snowpack properties)
Base.@kwdef struct ConstantDensity{Tρsn} <: SnowDensityScheme
    ρsn::Tρsn = 250.0u"kg/m^3" # constant snow density
end

function snowdensity!(
    snow::Snowpack{<:ConstantDensity},
    mass::DynamicSnowMassBalance,
    state
)
    ρsn = snow.para.density.ρsn
    ρw = waterdensity(snow)
    state.ρsn .= ρsn
    state.por .= 1 - ρsn / ρw
    return nothing
end
