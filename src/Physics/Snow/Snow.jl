module Snow

using CryoGrid: Top, SubSurfaceProcess, BoundaryProcess, BoundaryStyle, Dirichlet
using CryoGrid.InputOutput: Forcing
using CryoGrid.Physics.HeatConduction
using CryoGrid.Numerics
using CryoGrid.Utils

import CryoGrid
import CryoGrid.Physics
import CryoGrid.Physics.HeatConduction

using Unitful

export Snowpack, SnowThermalProperties

"""
Base type for snow paramterization schemes.
"""
abstract type SnowParameterization <: CryoGrid.Parameterization end

"""
    Bulk{Tdsn} <: SnowParameterization

Simple, bulk snow scheme where snowpack is represented as a single grid cell with homogenous state.
"""
Base.@kwdef struct Bulk{Tdsn,Tthresh} <: SnowParameterization
    dsn::Tdsn = 0.0u"m" # snow depth
    thresh::Tthresh = 0.05u"m" # snow threshold
end

Base.@kwdef struct SnowThermalProperties{Tρ,Tk,Tc} <: HeatConduction.ThermalProperties
    ρsn::Tρ = 350.0u"kg/m^3" # bulk snow density
    ksn::Tk = 0.5u"W/m/K" # thermal conductivity of fresh snow (Westermann et al. 2011)
    csn::Tc = 7.5e6u"J/K/m^3" # heat capacity of snow (Westermann et al. 2011)
end

struct SnowFree <: CryoGrid.Callback end
CryoGrid.CallbackStyle(::Type{SnowFree}) = CryoGrid.Continuous()

"""
    Snowpack{TPara<:SnowParameterization,TProp,TSp} <: CryoGrid.SubSurface

Generic representation of a ground surface snow pack.
"""
Base.@kwdef struct Snowpack{TPara<:SnowParameterization,TProp,TSp} <: CryoGrid.SubSurface
    para::TPara = Bulk()
    prop::TProp = SnowThermalProperties()
    sp::TSp = nothing
end

include("snow_bulk.jl")

end