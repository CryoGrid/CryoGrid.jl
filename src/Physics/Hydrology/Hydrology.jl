module Hydrology

import CryoGrid
import ConstructionBase

using CryoGrid
using CryoGrid.Physics
using CryoGrid.Numerics
using CryoGrid.Numerics: flux!, divergence!, ∇
using CryoGrid.Utils

using IfElse
using ModelParameters
using Unitful

export WaterBalance

Base.@kwdef struct HydraulicProperties{Tconsts,Tkwsat,Trb,Trc}
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
    WaterBalance{TFlow<:WaterFlow,Tdt,Tsp,TProp} <: CryoGrid.SubSurfaceProcess

Represents subsurface water transport processes.
"""
struct WaterBalance{TFlow<:WaterFlow,Tdt,Tsp,TProp} <: CryoGrid.SubSurfaceProcess
    flow::TFlow # vertical flow scheme
    prop::TProp # hydraulic parameters/constants
    dtlim::Tdt # dtlim
    sp::Tsp # user-defined specialization
end
"""
    kwsat(::SubSurface, water::WaterBalance)

Hydraulic conductivity at saturation.
"""
kwsat(::SubSurface, water::WaterBalance) = water.prop.kw_sat
"""
    maxwater(::SubSurface, ::WaterBalance, state, i)

Returns the maximum volumetric water content (saturation point) for grid cell `i`. Defaults to `1`.
"""
@inline maxwater(::SubSurface, ::WaterBalance, state, i) = one(eltype(state.sat))
"""
    watervariables(::WaterBalance)

Diagnostic variables shared by all implementations of WaterBalance.
"""
watervariables(::WaterBalance) = (
    Diagnostic(:θwi, OnGrid(Cells), domain=0..1), # total volumetric water+ice content
    Diagnostic(:θw, OnGrid(Cells), domain=0..1), # unfrozen/liquid volumetric water content
    Diagnostic(:θsat, OnGrid(Cells), domain=0..1), # maximum volumetric water content (saturation point)
    Diagnostic(:∂θwi∂t, OnGrid(Cells)), # divergence of total water content
    Diagnostic(:jw, OnGrid(Edges), u"m/s"), # water fluxes over grid cell boundaries
    Diagnostic(:kw, OnGrid(Edges), u"m/s", domain=0..Inf), # hydraulic conductivity (edges)
    Diagnostic(:kwc, OnGrid(Cells), u"m/s", domain=0..Inf), # hydraulic conductivity (cells)
)
# Helper methods
"""
    reductionfactor(water::WaterBalance, x)

Flux reduction factor for near-saturated conditions:
```math
r(x) = 1 - 1/(1+exp(-β(x-c)))
```
where β is a smoothness parameter and c is the "center" or shift parameter.
"""
reductionfactor(water::WaterBalance, x) = 1 - 1/(1+exp(-(x - water.prop.r_c)*water.prop.r_β))
"""
    resetfluxes!(::SubSurface, water::WaterBalance, state)

Resets flux terms (`jw` and `∂θwi∂t`) for `WaterBalance`.
"""
@inline function resetfluxes!(::SubSurface, water::WaterBalance, state)
    state.jw .= zero(eltype(state.jw))
    state.∂θwi∂t .= zero(eltype(state.∂θwi∂t))
end
function balancefluxes!(::SubSurface, water::WaterBalance, state)
    N = length(state.kw)
    state.jw[1] = min(max(state.jw[1], -state.θw[1]), state.θsat[1] - state.θwi[1])
    @inbounds for i in 2:N-1
        let θw_up = state.θw[i-1],
            θw_lo = state.θw[i],
            θwi_up = state.θwi[i-1],
            θwi_lo = state.θwi[i],
            θsat_up = state.θsat[i-1],
            θsat_lo = state.θsat[i],
            jw = state.jw[i];
            # limit flux based on
            # i) available water in cell above and
            # ii) free pore space in cell below
            max_flux_up = max(jw, θwi_up - θsat_up, -θw_lo) # upward flux is negative
            min_flux_down = min(jw, θsat_lo - θwi_lo, θw_up) # downward flux is positive
            # reduction factors
            r₁ = reductionfactor(water, state.sat[i-1])
            r₂ = reductionfactor(water, state.sat[i])
            state.jw[i] = r₁*max_flux_up*(jw < zero(jw)) + r₂*min_flux_down*(jw >= zero(jw))
        end
    end
    state.jw[end] = min(max(state.jw[end], state.θwi[end] - state.θsat[end]), state.θw[end])
end

export Rainfall, ConstantInfiltration
include("water_bc.jl")
export BucketScheme
include("water_bucket.jl")

# Constructors for WaterBalance
default_dtlim(::BucketScheme) = Physics.MaxDelta(0.1)
default_dtlim(::WaterFlow) = Physics.MaxDelta(Inf)
WaterBalance(flow::WaterFlow = BucketScheme(); prop = HydraulicProperties(), dtlim = default_dtlim(flow), sp = nothing) = WaterBalance(flow, prop, dtlim, sp)

end
