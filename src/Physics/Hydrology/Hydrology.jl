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
    WaterBalanceFormulation

Base type for different formulations of `WaterBalance`.
"""
abstract type WaterBalanceFormulation end
"""
    WaterBalance{TForm<:WaterBalanceFormulation,Tdt,TProp} <: CryoGrid.SubSurfaceProcess

Represents subsurface water transport processes.
"""
struct WaterBalance{TForm<:WaterBalanceFormulation,Tdt,TProp} <: CryoGrid.SubSurfaceProcess
    form::TForm
    prop::TProp
    dtlim::Tdt
end

export BucketScheme
include("water_bucket.jl")

# Constructors for WaterBalance
default_dtlim(::BucketScheme) = Physics.MaxDelta(0.1)
default_dtlim(::WaterBalanceFormulation) = Physics.MaxDelta(Inf)
WaterBalance(form::WaterBalanceFormulation = BucketScheme(); prop = Physics.Constants(), dtlim = default_dtlim(form)) = WaterBalance(form, prop, dtlim)

end
