abstract type SnowDensityScheme end

# constant density (using Snowpack properties)
Base.@kwdef struct ConstantDensity{Tρsn} <: SnowDensityScheme
    ρsn::Tρsn = 250.0u"kg/m^3" # constant snow density
end
