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
Base.@kwdef struct Bulk{Tdsn,Tρsn,Tthresh} <: SnowParameterization
    dsn::Tdsn = 0.0u"m" # snow depth
    ρsn::Tρsn = 250.0u"kg/m^3" # bulk snow density
    thresh::Tthresh = 0.05u"m" # snow threshold
end

struct SnowFree <: CryoGrid.Callback end
CryoGrid.CallbackStyle(::Type{SnowFree}) = CryoGrid.Continuous()

"""
    Snowpack{TPara<:SnowParameterization,TProp,TSp} <: CryoGrid.SubSurface

Generic representation of a ground surface snow pack.
"""
Base.@kwdef struct Snowpack{TPara<:SnowParameterization,TProp,TSp} <: CryoGrid.SubSurface
    para::TPara = Bulk()
    sp::TSp = nothing
end

include("snow_bulk.jl")

end