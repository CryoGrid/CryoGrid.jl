module Snow

using CryoGrid
using CryoGrid: ContinuousEvent, Increasing, Decreasing # for events/callbacks
using CryoGrid.InputOutput: Forcing
using CryoGrid.Physics.HeatConduction
using CryoGrid.Numerics
using CryoGrid.Utils

import CryoGrid
import CryoGrid.Physics
import CryoGrid.Physics.HeatConduction

using IfElse
using ModelParameters
using Unitful

export Snowpack, SnowProperties

SnowProperties(
    consts=Physics.Constants();
    ρw = consts.ρw,
    ρsn_new = Param(250.0, units=u"kg/m^3"),
    ρsn_old = Param(500.0, units=u"kg/m^3"),
) = (; ρw, ρsn_new, ρsn_old)

"""
    SnowpackParameterization

Base type for snowpack paramterization schemes.
"""
abstract type SnowpackParameterization <: CryoGrid.Parameterization end

"""
    Bulk{Tthresh} <: SnowpackParameterization

Simple, bulk ("single layer") snow scheme where snowpack is represented as a single grid cell with homogenous state.
"""
Base.@kwdef struct Bulk{Tthresh} <: SnowpackParameterization
    thresh::Tthresh = 0.02u"m" # snow threshold
end

"""
    Snowpack{Tpara<:SnowpackParameterization,Tprop,Tsp} <: CryoGrid.SubSurface

Generic representation of a ground surface snow pack.
"""
Base.@kwdef struct Snowpack{Tpara<:SnowpackParameterization,Tprop,Tsp} <: CryoGrid.SubSurface
    para::Tpara = Bulk()
    prop::Tprop = SnowProperties()
    sp::Tsp = nothing
end

abstract type SnowAblationScheme end
Base.@kwdef struct DegreeDayMelt{Tfactor,Tmax} <: SnowAblationScheme
    factor::Tfactor = 5.0u"mm/K/d"
    max_unfrozen::Tmax = 0.5
end

abstract type SnowAccumulationScheme end
Base.@kwdef struct LinearAccumulation{S} <: SnowAccumulationScheme
    rate_scale::S = Param(1.0, bounds=(0,Inf)) # scaling factor for snowfall rate
end

abstract type SnowDensityScheme end
# constant density (using Snowpack properties)
struct ConstantDensity end

abstract type SnowMassParameterization end
Base.@kwdef struct Prescribed{Tswe,Tρsn} <: SnowMassParameterization
    swe::Tswe = 0.0u"m" # depth snow water equivalent [m]
    ρsn::Tρsn = 250.0u"kg/m^3" # snow density [kg m^-3]
end
Base.@kwdef struct Dynamic{TAcc,TAbl,TDen} <: SnowMassParameterization
    accumulation::TAcc = LinearAccumulation()
    ablation::TAbl = DegreeDayMelt()
    density::TDen = ConstantDensity()
end

Base.@kwdef struct SnowMassBalance{Tpara} <: CryoGrid.SubSurfaceProcess
    para::Tpara = Dynamic()
end
accumulation(snow::SnowMassBalance{<:Dynamic}) = snow.para.accumulation
ablation(snow::SnowMassBalance{<:Dynamic}) = snow.para.ablation
density(snow::SnowMassBalance{<:Dynamic}) = snow.para.density

# convenience type aliases
"""
    PrescribedSnowMassBalance{Tswe,Tρsn} = SnowMassBalance{Prescribed{Tswe,Tρsn}} where {Tswe,Tρsn}
"""
const PrescribedSnowMassBalance{Tswe,Tρsn} = SnowMassBalance{Prescribed{Tswe,Tρsn}} where {Tswe,Tρsn}
"""
    DynamicSnowMassBalance{TAcc,TAbl,TDen} = SnowMassBalance{Dynamic{TAcc,TAbl,TDen}} where {TAcc,TAbl,TDen}
"""
const DynamicSnowMassBalance{TAcc,TAbl,TDen} = SnowMassBalance{Dynamic{TAcc,TAbl,TDen}} where {TAcc,TAbl,TDen}

snowvariables(::Snowpack, ::SnowMassBalance) = (
    Diagnostic(:dsn, Scalar, u"m"),
    Diagnostic(:T_ub, Scalar, u"°C"),
)

swe(::Snowpack, ::SnowMassBalance, state) = state.swe
swe(::Snowpack, smb::SnowMassBalance{<:Prescribed}, state) = smb.para.swe
swe(::Snowpack, smb::SnowMassBalance{<:Prescribed{<:Forcing}}, state) = smb.para.swe(state.t)
snowdensity(::Snowpack, ::SnowMassBalance, state) = state.ρsn
snowdensity(::Snowpack, smb::SnowMassBalance{<:Prescribed}, state) = smb.para.ρsn
snowdensity(::Snowpack, smb::SnowMassBalance{<:Prescribed{Tswe,<:Forcing}}, state) where {Tswe} = smb.para.ρsn(state.t)

# Boundary conditions

struct Snowfall{Tsn<:Forcing} <: BoundaryProcess{SnowMassBalance}
    snowfall::Tsn
end
CryoGrid.BoundaryStyle(::Snowfall) = CryoGrid.Neumann()
@inline boundaryvalue(bc::Snowfall, ::Top, ::SnowMassBalance, ::Snowpack, s1, s2) = bc.snowfall(s1.t)

# Implementations

# for prescribed snow depth/density, the mass balance is given so we do not need to do anything here
CryoGrid.prognosticstep!(::Snowpack, ::SnowMassBalance{<:Prescribed}, ssnow) = nothing

include("snow_bulk.jl")

end