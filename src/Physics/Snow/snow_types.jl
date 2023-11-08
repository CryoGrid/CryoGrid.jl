"""
    SnowpackParameterization

Base type for snowpack paramterization schemes.
"""
abstract type SnowpackParameterization <: CryoGrid.Parameterization end

"""
    SnowMassBalance{TAcc,TAbl} <: CryoGrid.SubSurfaceProcess

Subsurface process for snow layers governing how snow is accumulated and ablated.
"""
Base.@kwdef struct SnowMassBalance{TAcc,TAbl,TAux} <: CryoGrid.SubSurfaceProcess
    accumulation::TAcc = LinearAccumulation()
    ablation::TAbl = DegreeDayMelt()
    aux::TAux = nothing
end

"""
    SnowBC

Type alias for any `BoundaryProcess` compatible with `SnowMassBalance`.
"""
const SnowBC = BoundaryProcess{T} where {SnowMassBalance<:T<:SubSurfaceProcess}

"""
    SnowAblationScheme

Base type for different snow ablation (i.e. melting or redistribution) schemes.
"""
abstract type SnowAblationScheme end

"""
    SnowAccumulationScheme

Base type for different snow accumulation schemes.
"""
abstract type SnowAccumulationScheme end

"""
    SnowDensityScheme

Base type for different snow density schemes.
"""
abstract type SnowDensityScheme end

"""
    Snowpack{Tpara<:SnowpackParameterization,Tmass<:SnowMassBalance,Theat<:HeatBalance,Twater<:WaterBalance,Taux} <: CryoGrid.SubSurface

Generic representation of a snowpack "subsurface" layer.
"""
Base.@kwdef struct Snowpack{Tpara<:SnowpackParameterization,Tmass<:SnowMassBalance,Theat<:HeatBalance,Twater<:WaterBalance,Taux} <: CryoGrid.SubSurface
    para::Tpara = Bulk()
    mass::Tmass = SnowMassBalance()
    heat::Theat = HeatBalance()
    water::Twater = WaterBalance()
    aux::Taux = nothing
end

# Processes type aliases
const CoupledSnowWaterHeat{Tmass,Twater,Theat} = Coupled(SnowMassBalance, WaterBalance, HeatBalance)

"""
    Snowpack(para::SnowpackParameterization; kwargs...)

Convenience constructor that accepts the parameterization as a positional argument.
"""
Snowpack(para::SnowpackParameterization; kwargs...) = Snowpack(; para, kwargs...)

"""
    Bulk{Tden,Tthresh,Theat,Twater} <: SnowpackParameterization

Simple, bulk ("single layer") snow scheme where snowpack is represented as a single grid cell with homogenous state.
"""
Base.@kwdef struct Bulk{Tden<:SnowDensityScheme,Tthresh,Theat,Twater} <: SnowpackParameterization
    thresh::Tthresh = 0.02u"m" # snow threshold
    density::Tden = ConstantDensity() # snow density
    heat::Theat = ThermalProperties() # thermal properties
    water::Twater = HydraulicProperties(kw_sat=1e-4) # hydraulic properties
end

"""
    BulkSnowpack = Snowpack{<:Bulk}

Type alias for Snowpack with `Bulk` parameterization.
"""
const BulkSnowpack{TD} = Snowpack{<:Bulk{TD}} where {TD}
