Utils.@properties SnowMassProperties(
    ρw = CryoGrid.Constants.ρw,
    ρsn_new = 250.0u"kg/m^3",
    ρsn_old = 500.0u"kg/m^3",
)

"""
    SnowAblationScheme

Base type for different snow ablation (i.e. melting or redistribution) schemes.
"""
abstract type SnowAblationScheme end

Base.@kwdef struct DegreeDayMelt{Tfactor,Tmax} <: SnowAblationScheme
    factor::Tfactor = 5.0u"mm/K/d"
    max_unfrozen::Tmax = 0.5
end

abstract type SnowAccumulationScheme end

Base.@kwdef struct LinearAccumulation{S} <: SnowAccumulationScheme
    rate_scale::S = 1.0 # scaling factor for snowfall rate
end

abstract type SnowDensityScheme end

# constant density (using Snowpack properties)
struct ConstantDensity end

abstract type SnowMassParameterization end

Base.@kwdef struct PrescribedSnow{Tswe,Tρsn} <: SnowMassParameterization
    swe::Tswe = 0.0u"m" # depth snow water equivalent [m]
    ρsn::Tρsn = 250.0u"kg/m^3" # snow density [kg m^-3]
end

Base.@kwdef struct DynamicSnow{TAcc,TAbl,TDen} <: SnowMassParameterization
    accumulation::TAcc = LinearAccumulation()
    ablation::TAbl = DegreeDayMelt()
    density::TDen = ConstantDensity()
end

Base.@kwdef struct SnowMassBalance{Tpara} <: CryoGrid.SubSurfaceProcess
    para::Tpara = DynamicSnow()
end

# convenience type aliases
"""
    PrescribedSnowMassBalance{Tswe,Tρsn} = SnowMassBalance{PrescribedSnow{Tswe,Tρsn}} where {Tswe,Tρsn}
"""
const PrescribedSnowMassBalance{Tswe,Tρsn} = SnowMassBalance{PrescribedSnow{Tswe,Tρsn}} where {Tswe,Tρsn}
"""
    DynamicSnowMassBalance{TAcc,TAbl,TDen} = SnowMassBalance{DynamicSnow{TAcc,TAbl,TDen}} where {TAcc,TAbl,TDen}
"""
const DynamicSnowMassBalance{TAcc,TAbl,TDen} = SnowMassBalance{DynamicSnow{TAcc,TAbl,TDen}} where {TAcc,TAbl,TDen}

const SnowBC = BoundaryProcess{T} where {SnowMassBalance<:T<:SubSurfaceProcess}
