Base.@kwdef struct BucketScheme{Tfc,Tkwsat} <: WaterFluxesImpl
    θfc::Tfc = Param(0.5, domain=0..1)
    kw_sat::Tkwsat = Param(1e-5, domain=0..Inf, units=u"m/s")
end
fieldcapacity(::SubSurface, water::WaterFluxes{<:BucketScheme}) = water.impl.θfc
hydraulicconductivity(::SubSurface, water::WaterFluxes{<:BucketScheme}) = water.impl.kw_sat
# methods for computing diagnostic quantities
@inline function waterice!(sub::SubSurface, water::WaterFluxes, state)
    por = porosity(sub, state)
    @. state.θwi = state.sat*por
end
@inline function liquidwater!(::SubSurface, water::WaterFluxes, state)
    @. state.θw = state.θwi
end
@inline function hydraulicconductivity!(sub::SubSurface, water::WaterFluxes{<:BucketScheme}, state)
    kw_sat = hydraulicconductivity(sub, water)
    por = porosity(sub, state)
    @. state.kwc = kw_sat*state.θw/por
    Numerics.harmonicmean!(@view(state.kw[2:end-1]), state.kwc, Δ(state.grids.kw))
    state.kw[1] = state.kwc[1] 
    state.kw[end] = state.kwc[end]
end
@inline function resetfluxes!(::SubSurface, water::WaterFluxes, state)
    @. state.jw .= zero(eltype(state.jw))
end
# Helper methods
"""
    advectiveflux(θw_up, θwi_lo, θp_lo, θfc, kw)

Computes the advective downward water flux given liquid water content in the upper grid cell (θw_up),
total water/ice content in the lower grid cell (θwi_lo), total porespace in the lower grid cell (θp_lo),
the field capacity (θfc), and hydraulic conductivity at the boundary (kw).
"""
@inline function advectiveflux(θw_up, θwi_lo, θp_lo, θfc, kw)
    # downward advective flux due to gravity;
    # the grid in CryoGrid.jl decreases downward, so no negative sign is necessary
    jw = kw*(θw_up >= θfc)
    # limit flux based on i) available water in cell above and ii) free pore space in cell below
    return min(min(jw, max(θw_up, zero(jw))), θp_lo - θwi_lo)
end
# CryoGrid methods
CryoGrid.variables(::WaterFluxes{<:BucketScheme}) = (
    Prognostic(:sat, OnGrid(Cells), domain=0..1), # autmoatically generates dsat_wi and jsat_wi variables
    Diagnostic(:θwi, OnGrid(Cells), domain=0..1),
    Diagnostic(:θw, OnGrid(Cells), domain=0..1),
    Diagnostic(:jw, OnGrid(Edges), u"m/s"),
    Diagnostic(:kw, OnGrid(Edges), u"m/s", domain=0..Inf),
    Diagnostic(:kwc, OnGrid(Cells), u"m/s", domain=0..Inf),
)
CryoGrid.initialcondition!(sub::SubSurface, water::WaterFluxes, state) = CryoGrid.diagnosticstep!(sub, water, state)
function CryoGrid.diagnosticstep!(sub::SubSurface, water::WaterFluxes, state)
    resetfluxes!(sub, water, state)
    waterice!(sub, water, state)
    liquidwater!(sub, water, state)
    hydraulicconductivity!(sub, water, state)
end
function CryoGrid.prognosticstep!(sub::SubSurface, water::WaterFluxes{<:BucketScheme}, state)
    N = length(state.kw) # number of grid edges (including boundaries)
    # loop over grid
    @inbounds for i in 2:N-1 # note that the index is over grid *edges*
        let θwᵢ₋₁ = state.θw[i-1], # cell above edge i
            θwiᵢ = state.θwi[i], # cell below edge i
            θfc = fieldcapacity(sub, water),
            porᵢ₋₁ = porosity(sub, state, i-1),
            porᵢ = porosity(sub, state, i),
            kw = state.kw[i],
            δ = CryoGrid.thickness(sub, state, i);
            # compute fluxes over inner grid cell faces
            state.jw[i] = advectiveflux(θwᵢ₋₁, θwiᵢ, porᵢ, θfc, kw)
            # saturation flux for upper cell: ∂s∂t∂t = ∂θwi∂t/θp
            state.dsat[i-1] = -(state.jw[i] - state.jw[i-1]) / (δ * porᵢ₋₁)
        end
    end
    δₙ₋₁ = CryoGrid.thickness(sub, state, N-1)
    porₙ₋₁ = porosity(sub, state, N-1)
    state.dsat[end] = -(state.jw[end] - state.jw[end-1]) / (δₙ₋₁*porₙ₋₁)
end
function CryoGrid.interact!(sub1::SubSurface, water1::WaterFluxes{<:BucketScheme}, sub2::SubSurface, water2::WaterFluxes{<:BucketScheme}, state1, state2)
    θw₁ = state1.θw[end]
    θwi₂ = state2.θwi[1]
    θfc = fieldcapacity(sub1, water1) # take field capacity from upper layer where water would drain from
    por₂ = porosity(sub2, state2, 1)
    kwc₁ = state1.kwc[end]
    kwc₂ = state2.kwc[1]
    δ₁ = CryoGrid.thickness(sub1, state1, last)
    δ₂ = CryoGrid.thickness(sub2, state2, first)
    kw = state1.kw[end] = state2.kw[1] = Numerics.harmonicmean(kwc₁, kwc₂, δ₁, δ₂)
    # setting both jw[end] on the upper layer and jw[1] on the lower layer is redundant since they refer to the same
    # element of the same underlying state array, but it's nice for clarity
    state1.jw[end] = state2.jw[1] = advectiveflux(θw₁, θwi₂, por₂, θfc, kw)
    return nothing
end
function CryoGrid.timestep(::SubSurface, ::WaterFluxes{<:BucketScheme}, state)
    dtmax = Inf
    @inbounds for i in 1:length(state.sat)
        # solve for dt in:
        # sat + dt*∂sat∂t = 1 if ∂sat∂t > 0
        # sat + dt*∂sat∂t = 0 if ∂sat∂t <= 0
        # sets the maximum timestep to the dt which would saturate the grid cell;
        # will be Inf when ∂sat∂t = 0
        dt = IfElse.ifelse(
            sign(state.dsat[i]) > zero(eltype(state.dsat)),
            (one(eltype(state.sat)) - state.sat[i]) / state.dsat[i],
            -state.sat[i] / state.dsat[i]
        )
        dt = isfinite(dt) && dt > zero(dt) ? dt : Inf # make sure it's +Inf
        dtmax = min(dtmax, dt)
    end
    return dtmax
end
