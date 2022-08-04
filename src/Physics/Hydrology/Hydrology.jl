module Hydrology

import CryoGrid

using CryoGrid
using CryoGrid.Physics
using CryoGrid.Numerics
using CryoGrid.Numerics: flux!, divergence!

using IfElse
using ModelParameters
using Unitful

export WaterFluxes

abstract type WaterFluxesImpl end

struct WaterFluxes{TImpl<:WaterFluxesImpl,TProp} <: CryoGrid.SubSurfaceProcess
    impl::TImpl
    prop::TProp
end

export BucketScheme
include("water_bucket.jl")

# Constructors for WaterFluxes
WaterFluxes(impl::WaterFluxesImpl = BucketScheme(); prop = Physics.Constants()) = WaterFluxes(impl, prop)

end
