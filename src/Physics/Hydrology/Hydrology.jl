module Hydrology

import CryoGrid
import ConstructionBase

using CryoGrid
using CryoGrid.Physics
using CryoGrid.Numerics
using CryoGrid.Numerics: flux!, divergence!, ∇

using IfElse
using ModelParameters
using Unitful

export WaterBalance

"""
    WaterFlow

Base type for different formulations of water flow in `WaterBalance`.
"""
abstract type WaterFlow end
"""
    WaterBalance{TFlow<:WaterFlow,Tdt,TProp} <: CryoGrid.SubSurfaceProcess

Represents subsurface water transport processes.
"""
struct WaterBalance{TFlow<:WaterFlow,Tdt,TProp} <: CryoGrid.SubSurfaceProcess
    flow::TFlow
    prop::TProp
    dtlim::Tdt
end

"""
    watervariables(::WaterBalance)

Diagnostic variables shared by all implementations of WaterBalance.
"""
watervariables(::WaterBalance) = (
    Diagnostic(:θwi, OnGrid(Cells), domain=0..1),
    Diagnostic(:θw, OnGrid(Cells), domain=0..1),
    Diagnostic(:θsat, OnGrid(Cells), domain=0..1), # maximum volumetric water content (saturation point)
    Diagnostic(:dθwidt, OnGrid(Cells)),
    Diagnostic(:jw, OnGrid(Edges), u"m/s"),
    Diagnostic(:kw, OnGrid(Edges), u"m/s", domain=0..Inf),
    Diagnostic(:kwc, OnGrid(Cells), u"m/s", domain=0..Inf),
)

export BucketScheme
include("water_bucket.jl")

# Constructors for WaterBalance
default_dtlim(::BucketScheme) = Physics.MaxDelta(0.1)
default_dtlim(::WaterFlow) = Physics.MaxDelta(Inf)
WaterBalance(flow::WaterFlow = BucketScheme(); prop = Physics.Constants(), dtlim = default_dtlim(flow)) = WaterBalance(flow, prop, dtlim)

end
