Base.@kwdef struct BucketScheme{Tfc} <: WaterFlow
    fieldcap::Tfc = Param(0.2, domain=0..1)
end
fieldcapacity(::SubSurface, water::WaterBalance{<:BucketScheme}) = water.flow.fieldcap
# methods for computing diagnostic quantities
@inline function resetfluxes!(::SubSurface, water::WaterBalance, state)
    state.jw .= zero(eltype(state.jw))
    state.dθwidt .= zero(eltype(state.dθwidt))
end
@inline function waterice!(sub::SubSurface, water::WaterBalance{<:BucketScheme}, state)
    @inbounds for i in 1:length(state.sat)
        state.θsat[i] = saturation(sub, water, state, i)
        state.θwi[i] = state.sat[i]*state.θsat[i]
    end
end
@inline function liquidwater!(::SubSurface, water::WaterBalance, state)
    @. state.θw = state.θwi
end
@inline function hydraulicconductivity!(sub::SubSurface, water::WaterBalance{<:BucketScheme}, state)
    kw_sat = kwsat(sub, water)
    Δkw = Δ(state.grids.kw)
    @. state.kwc = kw_sat*state.θw / state.θsat
    state.kw[1] = state.kwc[1]
    state.kw[end] = state.kwc[end]
    Numerics.harmonicmean!(@view(state.kw[2:end-1]), state.kwc, Δkw)
end
# Helper methods
"""
    reductionfactor(::WaterBalance, x, c, β)

Flux reduction factor for near-saturated conditions:
```math
r(x) = 1 - 1/(1+exp(-β(x-c)))
```
where β is a smoothness parameter and c is the "center" or shift parameter.
"""
reductionfactor(::WaterBalance, x, c, β) = 1 - 1/(1+exp(-(x - c)β))
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
function balancefluxes!(sub::SubSurface, water::WaterBalance, state)
    N = length(state.kw)
    θmin = fieldcapacity(sub, water)
    state.jw[1] = min(max(state.jw[1], θmin - state.θw[1]), state.θsat[1] - state.θwi[1])
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
            max_flux_up = max(jw, θwi_up - θsat_up, θmin - θw_lo) # upward flux is negative
            min_flux_down = min(jw, θsat_lo - θwi_lo, θw_up - θmin) # downward flux is positive
            # reduction factors
            r₁ = reductionfactor(water, state.sat[i-1], 0.96, 1e3)
            r₂ = reductionfactor(water, state.sat[i], 0.96, 1e3)
            state.jw[i] = r₁*max_flux_up*(jw < zero(jw)) + r₂*min_flux_down*(jw >= zero(jw))
        end
    end
    state.jw[end] = min(max(state.jw[end], state.θwi[end] - state.θsat[end]), state.θw[end] - θmin)
end
# CryoGrid methods
CryoGrid.variables(water::WaterBalance{<:BucketScheme}) = (
    Prognostic(:sat, OnGrid(Cells), domain=0..1), # autmoatically generates dsat
    watervariables(water)..., # get common variables
)
function CryoGrid.initialcondition!(sub::SubSurface, water::WaterBalance, state)
    CryoGrid.diagnosticstep!(sub, water, state)
end
function CryoGrid.diagnosticstep!(sub::SubSurface, water::WaterBalance, state)
    resetfluxes!(sub, water, state)
    waterice!(sub, water, state)
    liquidwater!(sub, water, state)
    hydraulicconductivity!(sub, water, state)
end
function CryoGrid.prognosticstep!(sub::SubSurface, water::WaterBalance{<:BucketScheme}, state)
    wateradvection!(sub, water, state)
    balancefluxes!(sub, water, state)
    Numerics.divergence!(state.dθwidt, state.jw, Δ(state.grids.jw))
    @. state.dsat = state.dθwidt / state.θsat
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
    r₁ = reductionfactor(water1, state1.sat[end], 0.96, 1e3)
    r₂ = reductionfactor(water2, state2.sat[1], 0.96, 1e3)
    # setting both jw[end] on the upper layer and jw[1] on the lower layer is redundant since they refer to the same
    # element of the same underlying state array, but it's nice for clarity
    state1.jw[end] = state2.jw[1] = jw*r₁*(jw < zero(jw)) + jw*r₂*(jw >= zero(jw))
    return nothing
end
function CryoGrid.timestep(::SubSurface, water::WaterBalance{TFlow,TEvt,<:Physics.MaxDelta}, state) where {TFlow,TEvt}
    dtmax = Inf
    @inbounds for i in 1:length(state.sat)
        # solve for dt in:
        # sat + dt*∂sat∂t = 1 if ∂sat∂t > 0
        # sat + dt*∂sat∂t = 0 if ∂sat∂t <= 0
        # sets the maximum timestep to the dt which would saturate the grid cell;
        # will be Inf when ∂sat∂t = 0
        dt = water.dtlim(state.dsat[i], state.sat[i], state.t, zero(state.t), one(state.t))
        dt = isfinite(dt) && dt > zero(dt) ? dt : Inf # make sure it's +Inf
        dtmax = min(dtmax, dt)
    end
    return dtmax
end
