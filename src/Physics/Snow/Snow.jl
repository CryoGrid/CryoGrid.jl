module Snow

using CryoGrid
using CryoGrid: ContinuousEvent, Increasing, Decreasing # for events/callbacks
using CryoGrid.InputOutput: Forcing
using CryoGrid.Physics
using CryoGrid.Physics.Heat
using CryoGrid.Physics.Hydrology
using CryoGrid.Numerics
using CryoGrid.Utils

import CryoGrid
import CryoGrid.InputOutput
import CryoGrid.Physics
import CryoGrid.Physics.Heat

using IfElse
using ModelParameters
using Unitful
using UnPack

export Snowpack, SnowProperties, SnowMassBalance, Snowfall

SnowThermalProperties = Heat.ThermalProperties

Utils.@properties SnowProperties(
    ρw = Physics.Constants.ρw,
    ρsn_new = 250.0u"kg/m^3",
    ρsn_old = 500.0u"kg/m^3",
    heat = SnowThermalProperties(),
)

"""
    SnowpackParameterization

Base type for snowpack paramterization schemes.
"""
abstract type SnowpackParameterization <: CryoGrid.Parameterization end

"""
    Snowpack{Tpara<:SnowpackParameterization,Tprop,Tsp} <: CryoGrid.SubSurface

Generic representation of a ground surface snow pack.
"""
Base.@kwdef struct Snowpack{Tpara<:SnowpackParameterization,Tprop,Tsp,Tproc} <: CryoGrid.SubSurface{Tproc}
    para::Tpara = Bulk()
    prop::Tprop = SnowProperties()
    sp::Tsp = nothing
    proc::Tproc
end
Snowpack(proc::Tproc; kwargs...) where {Tproc} = Snowpack(;proc, kwargs...)

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
accumulation(snow::SnowMassBalance{<:DynamicSnow}) = snow.para.accumulation
ablation(snow::SnowMassBalance{<:DynamicSnow}) = snow.para.ablation
density(snow::SnowMassBalance{<:DynamicSnow}) = snow.para.density

# convenience type aliases
"""
    PrescribedSnowMassBalance{Tswe,Tρsn} = SnowMassBalance{PrescribedSnow{Tswe,Tρsn}} where {Tswe,Tρsn}
"""
const PrescribedSnowMassBalance{Tswe,Tρsn} = SnowMassBalance{PrescribedSnow{Tswe,Tρsn}} where {Tswe,Tρsn}
"""
    DynamicSnowMassBalance{TAcc,TAbl,TDen} = SnowMassBalance{DynamicSnow{TAcc,TAbl,TDen}} where {TAcc,TAbl,TDen}
"""
const DynamicSnowMassBalance{TAcc,TAbl,TDen} = SnowMassBalance{DynamicSnow{TAcc,TAbl,TDen}} where {TAcc,TAbl,TDen}

CryoGrid.basevariables(::Snowpack, ::SnowMassBalance) = (
    Diagnostic(:dsn, Scalar, u"m", domain=0..Inf),
    Diagnostic(:T_ub, Scalar, u"°C"),
)

swe(::Snowpack, ::SnowMassBalance, state) = state.swe
swe(::Snowpack, smb::SnowMassBalance{<:PrescribedSnow}, state) = smb.para.swe
swe(::Snowpack, smb::SnowMassBalance{<:PrescribedSnow{<:Forcing{u"m"}}}, state) = smb.para.swe(state.t)
snowdensity(::Snowpack, ::SnowMassBalance, state) = state.ρsn
snowdensity(::Snowpack, smb::SnowMassBalance{<:PrescribedSnow}, state) = smb.para.ρsn
snowdensity(::Snowpack, smb::SnowMassBalance{<:PrescribedSnow{Tswe,<:Forcing{u"kg/m^3"}}}, state) where {Tswe} = smb.para.ρsn(state.t)

# Boundary conditions
struct Snowfall{Tsn<:Forcing{u"m/s"}} <: BoundaryProcess{SnowMassBalance}
    snowfall::Tsn
end
CryoGrid.BoundaryStyle(::Snowfall) = CryoGrid.Neumann()
@inline boundaryvalue(bc::Snowfall, ::Top, ::SnowMassBalance, ::Snowpack, s1, s2) = bc.snowfall(s1.t)

# Implementations

# for prescribed snow depth/density, the mass balance is given so we do not need to do anything here
CryoGrid.prognosticstep!(::Snowpack, ::SnowMassBalance{<:PrescribedSnow}, ssnow) = nothing

# thermal properties snowpack
Heat.thermalproperties(snow::Snowpack) = snow.prop.heat
# volumetric fractions for snowpack
@inline function Physics.volumetricfractions(::Snowpack, state, i)
    @inbounds let θwi = state.θwi[i],
        θw = state.θw[i],
        θa = 1.0 - θwi,
        θi = θwi - θw;
        return (θw, θi, θa)
    end
end

include("snow_bulk.jl")

end