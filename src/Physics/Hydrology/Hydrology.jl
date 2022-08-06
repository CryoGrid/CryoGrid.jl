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

export WaterFluxes

"""
    WaterFluxesImpl

Base type for different formulations of `WaterFluxes`.
"""
abstract type WaterFluxesImpl end
"""
    WaterFluxes{TImpl<:WaterFluxesImpl,Tdt,TProp} <: CryoGrid.SubSurfaceProcess

Represents subsurface water transport processes.
"""
struct WaterFluxes{TImpl<:WaterFluxesImpl,Tdt,TProp} <: CryoGrid.SubSurfaceProcess
    impl::TImpl
    prop::TProp
    dtlim::Tdt
end

struct SaturationLimiter <: Physics.StepLimiter end

export BucketScheme
include("water_bucket.jl")

# Constructors for WaterFluxes
_default_dtlim(::BucketScheme) = SaturationLimiter()
_default_dtlim(::WaterFluxesImpl) = nothing
WaterFluxes(impl::WaterFluxesImpl = BucketScheme(:sat); prop = Physics.Constants(), dtlim = _default_dtlim(impl)) = WaterFluxes(impl, prop, dtlim)

end
