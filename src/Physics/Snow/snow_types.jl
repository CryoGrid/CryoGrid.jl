"""
    SnowpackParameterization

Base type for snowpack paramterization schemes.
"""
abstract type SnowpackParameterization <: CryoGrid.Parameterization end

"""
    SnowMassBalance{Tpara<:SnowMassParameterization} <: CryoGrid.SubSurfaceProcess

Base type for subsurface processes representing the dynamic accumulation and ablation of snow cover.
"""
abstract type SnowMassBalance <: CryoGrid.SubSurfaceProcess end

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
    PrescribedSnowMassBalance{Tswe} <: SnowMassBalance

"Prescribed" snow mass balance, i.e. where the snow water equivalent is given as a constant or forcing.
"""
Base.@kwdef struct PrescribedSnowMassBalance{Tswe} <: SnowMassBalance
    swe::Tswe = 0.0u"m" # depth snow water equivalent [m]
end

"""
    DynamicSnowMassBalance{TAcc,TAbl} <: SnowMassBalance

Dynamic snow mass balance, i.e. where snow is accumulated and ablated according to dynamic physical processes.
"""
Base.@kwdef struct DynamicSnowMassBalance{TAcc,TAbl} <: SnowMassBalance
    accumulation::TAcc = LinearAccumulation()
    ablation::TAbl = DegreeDayMelt()
end

"""
    Snowpack{Tpara<:SnowpackParameterization,Tmass<:SnowMassBalance,Theat<:HeatBalance,Twater<:WaterBalance,Taux} <: CryoGrid.SubSurface

Generic representation of a snowpack "subsurface" layer.
"""
Base.@kwdef struct Snowpack{Tpara<:SnowpackParameterization,Tmass<:SnowMassBalance,Theat<:HeatBalance,Twater<:WaterBalance,Taux} <: CryoGrid.SubSurface
    para::Tpara = Bulk()
    mass::Tmass = DynamicSnowMassBalance()
    heat::Theat = HeatBalance()
    water::Twater = WaterBalance()
    aux::Taux = nothing
end

# Processes type aliases
const CoupledSnowWaterHeat{Tmass,Twater,Theat} = Coupled(SnowMassBalance, WaterBalance, HeatBalance)
