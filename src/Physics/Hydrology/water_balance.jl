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
    BucketScheme{Tfc} <: WaterFlow

"Bucket" water scheme for downward advective flow due to gravity.
"""
Base.@kwdef struct BucketScheme{Tfc} <: WaterFlow
    fieldcap::Tfc = Param(0.2, domain=0..1)
end
# Constructors for WaterBalance
default_dtlim(::BucketScheme) = Physics.MaxDelta(0.1)
default_dtlim(::WaterFlow) = Physics.MaxDelta(Inf)
WaterBalance(flow::WaterFlow = BucketScheme(); prop = HydraulicProperties(), dtlim = default_dtlim(flow), sp = nothing) = WaterBalance(flow, prop, dtlim, sp)
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
CryoGrid.variables(water::WaterBalance) = (
    CryoGrid.variables(water.flow)...,
    CryoGrid.variables(water.evt)...,
    Diagnostic(:jw, OnGrid(Edges), u"m/s"), # water fluxes over grid cell boundaries
    Diagnostic(:jwET, OnGrid(Edges), u"m/s"), # water fluxes due to evapotranspiration
    Diagnostic(:θwi, OnGrid(Cells), domain=0..1), # total volumetric water+ice content
    Diagnostic(:θw, OnGrid(Cells), domain=0..1), # unfrozen/liquid volumetric water content
    Diagnostic(:θsat, OnGrid(Cells), domain=0..1), # maximum volumetric water content (saturation point)
    Diagnostic(:∂θwi∂t, OnGrid(Cells)), # divergence of total water content
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
    state.jwET .= zero(eltype(state.jwET))
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
fieldcapacity(::SubSurface, water::WaterBalance{<:BucketScheme}) = water.flow.fieldcap
@inline function watercontent!(sub::SubSurface, water::WaterBalance{<:BucketScheme}, state)
    @inbounds for i in 1:length(state.sat)
        state.θsat[i] = maxwater(sub, water, state, i)
        state.θwi[i] = state.sat[i]*state.θsat[i]
    end
end
@inline function hydraulicconductivity!(sub::SubSurface, water::WaterBalance{<:BucketScheme}, state)
    kw_sat = kwsat(sub, water)
    Δkw = Δ(state.grids.kw)
    @. state.kwc = kw_sat*state.θw / state.θsat
    state.kw[1] = state.kwc[1]
    state.kw[end] = state.kwc[end]
    Numerics.harmonicmean!(@view(state.kw[2:end-1]), state.kwc, Δkw)
end
"""
    advectiveflux(θw_up, θwi_lo, θsat_lo, θmin, kw)

Computes the advective downward water flux given liquid water content in the upper grid cell (θw_up),
the minimum water content (θmin), and hydraulic conductivity at the boundary (kw).
"""
@inline function advectiveflux(θw_up, θmin, kw)
    # downward advective flux due to gravity;
    # the grid in CryoGrid.jl is positive downward, so no negative sign is necessary
    jw = kw*(θw_up > θmin)
    return jw
end
"""
    wateradvection!(sub::SubSurface, water::WaterBalance, state)

Computes the advective component of water fluxes due to gravity and stores the result in `state.jw`.
"""
function wateradvection!(sub::SubSurface, water::WaterBalance, state)
    N = length(state.kw) # number of grid edges (including boundaries)
    # loop over grid
    @inbounds for i in 2:N-1 # note that the index is over grid *edges*
        let θwᵢ₋₁ = state.θw[i-1], # cell above edge i
            θfc = fieldcapacity(sub, water),
            kw = state.kw[i];
            # compute fluxes over inner grid cell faces
            state.jw[i] = advectiveflux(θwᵢ₋₁, θfc, kw)
        end
    end
end
# CryoGrid methods
CryoGrid.variables(::BucketScheme) = (
    Prognostic(:sat, OnGrid(Cells), domain=0..1), # autmoatically generates ∂sat∂t
)
function CryoGrid.initialcondition!(sub::SubSurface, water::WaterBalance, state)
    CryoGrid.diagnosticstep!(sub, water, state)
end
function CryoGrid.diagnosticstep!(sub::SubSurface, water::WaterBalance, state)
    resetfluxes!(sub, water, state)
    watercontent!(sub, water, state)
    hydraulicconductivity!(sub, water, state)
end
function CryoGrid.prognosticstep!(sub::SubSurface, water::WaterBalance{<:BucketScheme}, state)
    wateradvection!(sub, water, state)
    balancefluxes!(sub, water, state)
    Numerics.divergence!(state.∂θwi∂t, state.jw, Δ(state.grids.jw))
    @. state.∂sat∂t = state.∂θwi∂t / state.θsat
end
function CryoGrid.interact!(sub1::SubSurface, water1::WaterBalance{<:BucketScheme}, sub2::SubSurface, ::WaterBalance{<:BucketScheme}, state1, state2)
    θw₁ = state1.θw[end]
    θfc = fieldcapacity(sub1, water1) # take field capacity from upper layer where water would drain from
    kwc₁ = state1.kwc[end]
    kwc₂ = state2.kwc[1]
    δ₁ = CryoGrid.thickness(sub1, state1, last)
    δ₂ = CryoGrid.thickness(sub2, state2, first)
    kw = state1.kw[end] = state2.kw[1] = Numerics.harmonicmean(kwc₁, kwc₂, δ₁, δ₂)
    jw = advectiveflux(θw₁, θfc, kw)
    # reduction factors
    r₁ = reductionfactor(water1, state1.sat[end])
    r₂ = reductionfactor(water2, state2.sat[1])
    # setting both jw[end] on the upper layer and jw[1] on the lower layer is redundant since they refer to the same
    # element of the same underlying state array, but it's nice for clarity
    state1.jw[end] = state2.jw[1] = jw*r₁*(jw < zero(jw)) + jw*r₂*(jw >= zero(jw))
    return nothing
end
function CryoGrid.timestep(::SubSurface, water::WaterBalance{TFlow,<:Physics.MaxDelta}, state) where {TFlow}
    dtmax = Inf
    @inbounds for i in 1:length(state.sat)
        # solve for dt in:
        # sat + dt*∂sat∂t = 1 if ∂sat∂t > 0
        # sat + dt*∂sat∂t = 0 if ∂sat∂t <= 0
        # sets the maximum timestep to the dt which would saturate the grid cell;
        # will be Inf when ∂sat∂t = 0
        dt = water.dtlim(state.∂sat∂t[i], state.sat[i], state.t, zero(state.t), one(state.t))
        dt = isfinite(dt) && dt > zero(dt) ? dt : Inf # make sure it's +Inf
        dtmax = min(dtmax, dt)
    end
    return dtmax
end
