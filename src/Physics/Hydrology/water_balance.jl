# Default minwater for bucket scheme (must be nonzero for ET)
minwater(::SubSurface, water::WaterBalance{<:BucketScheme}) = 0.01

function limit_upper_flux(water::WaterBalance, jw, θw, θwi, θsat, sat, Δz)
    # case (i): jw < 0 -> outflow -> limit based on available water
    # case (ii): jw > 0 -> inflow -> limit based on available space
    # case (iii): jw = 0 -> no flow, limit has no effect
    jw = min(max(jw, -θw*Δz), (θsat - θwi)*Δz)
    # influx reduction factor
    r = reductionfactor(water, sat)*(jw > 0) + 1.0*(jw <= 0)
    return r*jw
end

function limit_lower_flux(water::WaterBalance, jw, θw, θwi, θsat, sat, Δz)
    # case (i): jw < 0 -> inflow -> limit based on available space
    # case (ii): jw > 0 -> outflow -> limit based on available water
    # case (iii): jw = 0 -> no flow, limit has no effect
    jw = min(max(jw, (θwi - θsat)*Δz), θw*Δz)
    # influx reduction factor
    r = reductionfactor(water, sat)*(jw < 0) + 1.0*(jw >= 0)
    return r*jw
end

"""
    balanceflux([water_up::WaterBalance], water_lo::WaterBalance, jw, θw_up, θw_lo, θwi_up, θwi_lo, θsat_up, θsat_lo, sat_up, sat_lo, Δz_up, Δz_lo)

Rescales the water flux jw to satisfy mass conservation based on the given hydrological states for two adjacent
grid cells ('up' refers to upper and 'lo' refers to lower). Also applies nonlinear reduction factors for near-saturated
conditions.
"""
balanceflux(water::WaterBalance, args...) = balanceflux(water, water, args...)
function balanceflux(water_up::WaterBalance, water_lo::WaterBalance, jw, θw_up, θw_lo, θwi_up, θwi_lo, θsat_up, θsat_lo, sat_up, sat_lo, Δz_up, Δz_lo)
    # flux for upper boundary of lower grid cell
    jw_lo = limit_upper_flux(water_lo, jw, θw_lo, θwi_lo, θsat_lo, sat_lo, Δz_lo)
    # flux for lower boundary of upper grid cell
    jw_up = limit_lower_flux(water_up, jw, θw_up, θwi_up, θsat_up, sat_up, Δz_up)
    # take min if positive flux, max otherwise
    return max(jw, jw_lo, jw_up)*(jw <= 0) + min(jw, jw_lo, jw_up)*(jw > 0)
end

"""
    balancefluxes!(sub::SubSurface, water::WaterBalance, state)

Applies `balancefluxes!` to all internal cell faces.
"""
function balancefluxes!(sub::SubSurface, water::WaterBalance, state)
    N = length(state.kw)
    dt = state.dt
    @inbounds for i in 2:N-1
        let θw_up = state.θw[i-1],
            θw_lo = state.θw[i],
            θwi_up = state.θwi[i-1],
            θwi_lo = state.θwi[i],
            θsat_up = state.θsat[i-1],
            θsat_lo = state.θsat[i],
            sat_up = state.sat[i-1],
            sat_lo = state.sat[i],
            Δz_up = CryoGrid.thickness(sub, state, i-1),
            Δz_lo = CryoGrid.thickness(sub, state, i),
            jw = state.jw[i]*dt;
            state.jw[i] = balanceflux(water, jw, θw_up, θw_lo, θwi_up, θwi_lo, θsat_up, θsat_lo, sat_up, sat_lo, Δz_up, Δz_lo)
        end
    end
    return nothing
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
    @inbounds for i in 2:N-1 # note that the index is over grid edges
        let θwᵢ₋₁ = state.θw[i-1], # cell above edge i
            θfc = minwater(sub, water, state, i),
            kw = state.kw[i];
            # compute fluxes over inner grid cell faces
            state.jw[i] += advectiveflux(θwᵢ₋₁, θfc, kw)
        end
    end
end

"""
    waterdiffusion!(::SubSurface, ::WaterBalance, state)

Computes diffusive fluxes for water balance, if defined. Default implementation does nothing.
"""
waterdiffusion!(::SubSurface, ::WaterBalance, state) = nothing

"""
    waterprognostic!(::SubSurface, ::WaterBalance, state)

Computes the prognostic time derivative for the water balance, usually based on `∂θwi∂t`.
Implementation depends on which water flow scheme is being used.
"""
function waterprognostic!(::SubSurface, ::WaterBalance{<:BucketScheme}, state)
    @inbounds @. state.∂sat∂t = state.∂θwi∂t / state.θsat
end

# Helper methods
"""
    reductionfactor(water::WaterBalance, x)

Hydraulic conductivity reduction factor for near-saturated conditions:
```math
r(x) = (1-exp(-βx^2))*(1-exp(-β*(1-x)^2))
```
where β is a smoothness parameter and c is the "center" or shift parameter.
"""
reductionfactor(water::WaterBalance, x) = 1 - exp(-water.prop.r_β*(1-x)^2)

"""
    resetfluxes!(::SubSurface, water::WaterBalance, state)

Resets flux terms (`jw` and `∂θwi∂t`) for `WaterBalance`.
"""
@inline function resetfluxes!(::SubSurface, water::WaterBalance, state)
    state.jw .= zero(eltype(state.jw))
    state.jw_ET .= zero(eltype(state.jw_ET))
    state.∂θwi∂t .= zero(eltype(state.∂θwi∂t))
end

# CryoGrid methods
CryoGrid.variables(water::WaterBalance) = (
    CryoGrid.variables(water.flow)...,
    CryoGrid.variables(water.et)...,
    Diagnostic(:jw, OnGrid(Edges), u"m/s"), # water fluxes over grid cell boundaries
    Diagnostic(:jw_ET, OnGrid(Edges), u"m/s"), # water fluxes due to evapotranspiration
    Diagnostic(:θwi, OnGrid(Cells), domain=0..1), # total volumetric water+ice content
    Diagnostic(:θw, OnGrid(Cells), domain=0..1), # unfrozen/liquid volumetric water content
    Diagnostic(:θsat, OnGrid(Cells), domain=0..1), # maximum volumetric water content (saturation point)
    Diagnostic(:∂θwi∂t, OnGrid(Cells)), # divergence of total water content
    Diagnostic(:kw, OnGrid(Edges), u"m/s", domain=0..Inf), # hydraulic conductivity (edges)
    Diagnostic(:kwc, OnGrid(Cells), u"m/s", domain=0..Inf), # hydraulic conductivity (cells)
)
CryoGrid.variables(::BucketScheme) = (
    Prognostic(:sat, OnGrid(Cells), domain=0..1), # autmoatically generates ∂sat∂t
)

function CryoGrid.initialcondition!(sub::SubSurface, water::WaterBalance, state)
    CryoGrid.updatestate!(sub, water, state)
end

function CryoGrid.updatestate!(sub::SubSurface, water::WaterBalance, state)
    resetfluxes!(sub, water, state)
    watercontent!(sub, water, state)
    hydraulicconductivity!(sub, water, state)
    evapotranspiration!(sub, water, state)
end

function CryoGrid.computefluxes!(sub::SubSurface, water::WaterBalance, state)
    evapotranspirative_fluxes!(sub, water, state)
    wateradvection!(sub, water, state)
    balancefluxes!(sub, water, state)
    Numerics.divergence!(state.∂θwi∂t, state.jw, Δ(state.grid))
    waterprognostic!(sub, water, state)
end

function CryoGrid.interact!(sub1::SubSurface, water1::WaterBalance{<:BucketScheme}, sub2::SubSurface, water2::WaterBalance{<:BucketScheme}, state1, state2)
    θw₁ = state1.θw[end]
    θw₂ = state2.θw[1]
    θwi₁ = state1.θwi[end]
    θwi₂ = state2.θwi[1]
    θsat₁ = state1.θsat[end]
    θsat₂ = state2.θsat[1]
    sat₁ = state1.sat[end]
    sat₂ = state2.sat[1]
    Δz₁ = CryoGrid.thickness(sub1, state1, last)
    Δz₂ = CryoGrid.thickness(sub2, state2, first)
    # take minimum water content from upper layer where water would drain from
    θmin₁ = minwater(sub1, water1, state1, lastindex(state1.θw))
    kwc₁ = state1.kwc[end]
    kwc₂ = state2.kwc[1]
    kw = state1.kw[end] = state2.kw[1] = harmonicmean(kwc₁, kwc₂, Δz₁, Δz₂)
    jw = advectiveflux(θw₁, θmin₁, kw)*state1.dt
    # setting both jw[end] on the upper layer and jw[1] on the lower layer is redundant since they refer to the same
    # element of the same underlying state array, but it's nice for clarity
    state1.jw[end] = state2.jw[1] = balanceflux(water1, water2, jw, θw₁, θw₂, θwi₁, θwi₂, θsat₁, θsat₂, sat₁, sat₂, Δz₁, Δz₂)
    return nothing
end

function CryoGrid.timestep(
    sub::SubSurface,
    water::WaterBalance{TFlow,TET,<:CryoGrid.MaxDelta},
    state
) where {TFlow,TET}
    dtmax = Inf
    @inbounds for i in 1:length(state.sat)
        dt = water.dtlim(state.∂sat∂t[i], state.sat[i], state.t, zero(state.t), one(state.t))
        dt = isfinite(dt) && dt > zero(dt) ? dt : Inf # make sure it's +Inf
        dtmax = min(dtmax, dt)
    end
    return dtmax
end

# No flow case
resetfluxes!(::SubSurface, ::WaterBalance{NoFlow}, state) = nothing

hydraulicconductivity!(::SubSurface, ::WaterBalance{NoFlow}, state) = nothing

CryoGrid.variables(::NoFlow) = (
    Diagnostic(:sat, OnGrid(Cells), domain=0..1), # autmoatically generates ∂sat∂t
)

function CryoGrid.initialcondition!(sub::SubSurface, water::WaterBalance{NoFlow}, state)
    @inbounds for i in eachindex(state.sat)
        state.θwi[i] = state.sat[i]*maxwater(sub, water, state, i)
    end
end

CryoGrid.computefluxes!(::SubSurface, ::WaterBalance{NoFlow}, state) = nothing

CryoGrid.updatestate!(::SubSurface, ::WaterBalance{NoFlow}, state) = nothing

CryoGrid.timestep(::SubSurface, ::WaterBalance{NoFlow}, state) = Inf
