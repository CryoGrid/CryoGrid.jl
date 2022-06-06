module Hydrology

import CryoGrid

using CryoGrid
using CryoGrid.Physics
using CryoGrid.Numerics
using CryoGrid.Numerics: nonlineardiffusion!

using IfElse
using ModelParameters
using Unitful

abstract type WaterFlowParameterization end

Base.@kwdef struct WaterFlow{TPara<:WaterFlowParameterization,Tinit,TProp} <: CryoGrid.SubSurfaceProcess
    para::TPara = BucketScheme()
    prop::TProp = HydrologicalProperties()
    init::Tinit = nothing # optional initialization scheme
end

include("water_bucket.jl")

end
