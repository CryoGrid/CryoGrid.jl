"""
    kwsat(::SubSurface, ::WaterBalance)

Hydraulic conductivity at saturation.
"""
kwsat(sub::SubSurface, ::WaterBalance) = hydraulicproperties(sub).kw_sat
"""
    maxwater(::SubSurface, ::WaterBalance, state, i=nothing)

Returns the maximum volumetric water content (saturation point) for grid cell `i`. Defaults to `1`.
"""
maxwater(::SubSurface, ::WaterBalance, state, i=nothing) = one(eltype(state.sat))
"""
    minwater(::SubSurface, water::WaterBalance, i=nothing)

Returns the minimum volumetric water content (typically field capacity for simplified schemes) for grid cell `i`. Defaults to zero.
"""
minwater(::SubSurface, water::WaterBalance, i=nothing) = 0.0
minwater(::SubSurface, water::WaterBalance{<:BucketScheme}, i=nothing) = water.flow.fieldcap
function balancefluxes!(::SubSurface, water::WaterBalance, state)
    # top flux
    r_top = reductionfactor(water, state.sat[1])
    state.jw[1] = r_top*min(max(state.jw[1], -state.θw[1]), state.θsat[1] - state.θwi[1])
    # inner cell faces
    @inbounds for i in 2:length(state.kw)-1
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
            state.jw[i] = max_flux_up*(jw < zero(jw)) + min_flux_down*(jw >= zero(jw))
            # apply reduction factors
            r_up = reductionfactor(water, state.sat[i-1])
            r_lo = reductionfactor(water, state.sat[i])
            state.jw[i] *= r_up*r_lo
        end
    end
    # bottom flux
    r_bot = reductionfactor(water, state.sat[end])
    state.jw[end] = r_bot*min(max(state.jw[end], state.θwi[end] - state.θsat[end]), state.θw[end])
end
watercontent(::SubSurface, state) = state.θwi
watercontent(sub::SubSurface, state, i) = Utils.getscalar(watercontent(sub, state), i)
@inline function watercontent!(sub::SubSurface, water::WaterBalance, state)
    @inbounds for i in 1:length(state.sat)
        state.θsat[i] = maxwater(sub, water, state, i)
        state.θwi[i] = state.sat[i]*state.θsat[i]
    end
end
@inline function hydraulicconductivity!(sub::SubSurface, water::WaterBalance{<:BucketScheme}, state)
    kw_sat = kwsat(sub, water)
    Δkw = Δ(state.grids.kw)
    state.kwc[1] = kw_sat*state.θw[1] / state.θsat[1]
    state.kw[1] = state.kwc[1]*reductionfactor(water, state.sat[1])
    # inner edges
    @inbounds for i in 2:length(state.kwc)-1
        state.kwc[i] = kw_sat*state.θw[i] / state.θsat[i]
        state.kw[i] = Numerics.harmonicmean(state.kwc[i-1], state.kwc[i], Δkw[i-1], Δkw[i])
        state.kw[i] *= reductionfactor(water, state.sat[i-1])*reductionfactor(water, state.sat[i])
    end
    state.kwc[end] = kw_sat*state.θw[end] / state.θsat[end]
    state.kw[end] = state.kwc[end]*reductionfactor(water, state.sat[end])
end
function wateradvection!(sub::SubSurface, water::WaterBalance, state)
    N = length(state.kw) # number of grid edges (including boundaries)
    # loop over grid
    @inbounds for i in 2:N-1 # note that the index is over grid edges
        let θwᵢ₋₁ = state.θw[i-1], # cell above edge i
            θfc = minwater(sub, water),
            kw = state.kw[i];
            # compute fluxes over inner grid cell faces
            state.jw[i] += advectiveflux(θwᵢ₋₁, θfc, kw)
        end
    end
end
waterdiffusion!(::SubSurface, ::WaterBalance, state) = nothing
function waterprognostic!(::SubSurface, ::WaterBalance{<:BucketScheme}, state)
    @inbounds @. state.∂sat∂t = state.∂θwi∂t / state.θsat
end
# Helper methods
"""
    reductionfactor(water::WaterBalance, x)

Hydraulic conductivity reduction factor for near-dry and near-saturated conditions:
```math
r(x) = (1-exp(-βx^2))*(1-exp(-β*(1-x)^2))
```
where β is a smoothness parameter and c is the "center" or shift parameter.
"""
reductionfactor(water::WaterBalance, x) = (1-exp(-water.prop.r_β*x^2))*(1-exp(-water.prop.r_β*(1-x)^2))
"""
    resetfluxes!(::SubSurface, water::WaterBalance, state)

Resets flux terms (`jw` and `∂θwi∂t`) for `WaterBalance`.
"""
@inline function resetfluxes!(::SubSurface, water::WaterBalance, state)
    state.jw .= zero(eltype(state.jw))
    state.jwET .= zero(eltype(state.jwET))
    state.∂θwi∂t .= zero(eltype(state.∂θwi∂t))
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
# CryoGrid methods
CryoGrid.variables(water::WaterBalance) = (
    CryoGrid.variables(water.flow)...,
    CryoGrid.variables(water.et)...,
    Diagnostic(:jw, OnGrid(Edges), u"m/s"), # water fluxes over grid cell boundaries
    Diagnostic(:jwET, OnGrid(Edges), u"m/s"), # water fluxes due to evapotranspiration
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
    CryoGrid.diagnosticstep!(sub, water, state)
end
function CryoGrid.diagnosticstep!(sub::SubSurface, water::WaterBalance, state)
    resetfluxes!(sub, water, state)
    watercontent!(sub, water, state)
    hydraulicconductivity!(sub, water, state)
    evapotranspiration!(sub, water, state)
end
function CryoGrid.prognosticstep!(sub::SubSurface, water::WaterBalance, state)
    evapotranspirative_fluxes!(sub, water, state)
    wateradvection!(sub, water, state)
    balancefluxes!(sub, water, state)
    Numerics.divergence!(state.∂θwi∂t, state.jw, Δ(state.grids.jw))
    waterprognostic!(sub, water, state)
end
function CryoGrid.interact!(sub1::SubSurface, water1::WaterBalance{<:BucketScheme}, sub2::SubSurface, water2::WaterBalance{<:BucketScheme}, state1, state2)
    θw₁ = state1.θw[end]
    θfc = minwater(sub1, water1) # take field capacity from upper layer where water would drain from
    kwc₁ = state1.kwc[end]
    kwc₂ = state2.kwc[1]
    δ₁ = CryoGrid.thickness(sub1, state1, last)
    δ₂ = CryoGrid.thickness(sub2, state2, first)
    kw = state1.kw[end] = state2.kw[1] = harmonicmean(kwc₁, kwc₂, δ₁, δ₂)
    jw = advectiveflux(θw₁, θfc, kw)
    # setting both jw[end] on the upper layer and jw[1] on the lower layer is redundant since they refer to the same
    # element of the same underlying state array, but it's nice for clarity
    state1.jw[end] = state2.jw[1] = jw*(jw < zero(jw)) + jw*(jw >= zero(jw))
    return nothing
end
function CryoGrid.timestep(
    ::SubSurface,
    water::WaterBalance{TFlow,TET,<:Physics.MaxDelta},
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
Hydrology.resetfluxes!(::SubSurface, ::WaterBalance{NoFlow}, state) = nothing
Hydrology.hydraulicconductivity!(::SubSurface, ::WaterBalance{NoFlow}, state) = nothing
CryoGrid.variables(::NoFlow) = (
    Diagnostic(:sat, OnGrid(Cells), domain=0..1), # autmoatically generates ∂sat∂t
)
function CryoGrid.initialcondition!(sub::SubSurface, water::WaterBalance{NoFlow}, state)
    @inbounds for i in eachindex(state.sat)
        state.θwi[i] = state.sat[i]*maxwater(sub, water, state, i)
    end
end
CryoGrid.prognosticstep!(::SubSurface, ::WaterBalance{NoFlow}, state) = nothing
CryoGrid.diagnosticstep!(::SubSurface, ::WaterBalance{NoFlow}, state) = nothing
CryoGrid.timestep(::SubSurface, ::WaterBalance{NoFlow}, state) = Inf
