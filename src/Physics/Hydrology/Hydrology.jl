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

Base.@kwdef struct HydrualicProperties{Tconsts,Tkwsat,Trb,Trc}
    consts::Tconsts = Physics.Constants()
    kw_sat::Tkwsat = Param(1e-5, domain=0..Inf, units=u"m/s")
    r_β::Trb = 1e3 # reduction factor scale parameter
    r_c::Trc = 0.96325 # reduction factor shift parameter
end

"""
    WaterFlow

Base type for different formulations of water flow in `WaterBalance`.
"""
abstract type WaterFlow end
"""
    WaterBalance{TFlow<:WaterFlow,TEvt,Tdt,TProp} <: CryoGrid.SubSurfaceProcess

Represents subsurface water transport processes.
"""
struct WaterBalance{TFlow<:WaterFlow,TEvt,Tdt,TProp} <: CryoGrid.SubSurfaceProcess
    flow::TFlow
    evt::TEvt
    prop::TProp
    dtlim::Tdt
end
"""
    kwsat(::SubSurface, water::WaterBalance)

Hydraulic conductivity at saturation.
"""
kwsat(::SubSurface, water::WaterBalance) = water.prop.kw_sat
"""
    evapotranspiration(::SubSurface, water::WaterBalance)

Evapotranspiration implementation.
"""
evapotranspiration(::SubSurface, water::WaterBalance) = water.evt
"""
    saturation(sub::SubSurface, ::WaterBalance, state, i)

Returns the maximum volumetric water content (saturation point) for grid cell `i`. Defaults to `1.0`.
"""
@inline saturation(sub::SubSurface, ::WaterBalance, state, i) = 1.0

"""
    Evapotranspiration

Base type for various parameterizations of evapotranspiration.
"""
abstract type Evapotranspiration end
"""
    SurfaceEvaporation <: Evapotranspiration

Evaporation-only scheme which computes fluxes for the top-most grid cell when interacting with a boundary layer.
"""
struct SurfaceEvaporation <: Evapotranspiration end
"""
    DampedET

Corresponds to evapotranspiration scheme 3 described in section 2.2.2 of Westermann et al. (2022).
"""
Base.@kwdef struct DampedET{Tftr,Tfev,Tdtr,Tdev} <: Evapotranspiration
    f_tr::Tftr = Param(0.5)
    f_ev::Tfev = Param(0.5)
    d_tr::Tdtr = Param(1.0, units=u"m")
    d_ev::Tdev = Param(0.1, units=u"m")
end

# auxiliary flux variables
auxiliaryfluxes(::WaterBalance) = ()
auxiliaryfluxes(::WaterBalance{<:WaterFlow,<:Evapotranspiration}) = (
    Diagnostic(:jwET, OnGrid(Edges), u"m/s"),
)

"""
    watervariables(::WaterBalance)

Diagnostic variables shared by all implementations of WaterBalance.
"""
watervariables(water::WaterBalance) = (
    auxiliaryfluxes(water)...,
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
WaterBalance(flow::WaterFlow = BucketScheme(); prop = HydrualicProperties(), evt = nothing, dtlim = default_dtlim(flow)) = WaterBalance(flow, prop, evt, dtlim)

end
