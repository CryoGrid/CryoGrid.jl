"""
    SnowpackParameterization

Base type for snowpack paramterization schemes.
"""
abstract type SnowpackParameterization <: CryoGrid.Parameterization end

"""
    SnowMassBalance{TAcc,TAbl} <: CryoGrid.SubSurfaceProcess

Subsurface process for snow layers governing how snow is accumulated and ablated.
"""
Base.@kwdef struct SnowMassBalance{TAcc,TAbl,TAux,TDt} <: CryoGrid.SubSurfaceProcess
    accumulation::TAcc = LinearAccumulation()
    ablation::TAbl = DegreeDayMelt()
    dtlim::TDt = CryoGrid.MaxDelta(0.01)
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
    SnowThermalConductivity

Base type for snow thermal conductivity parameterizations.
"""
abstract type SnowThermalConductivity end

"""
    SnowThermalProperties{Tcond<:SnowThermalConductivity,Tprop}

Specifies the thermal properties of the snowpack.
"""
Base.@kwdef struct SnowThermalProperties{Tcond<:SnowThermalConductivity,Tprop}
    cond::Tcond = SturmQuadratic()
    prop::Tprop = ThermalProperties()
end

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
