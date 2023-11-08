"""
    watercontent!(::SubSurface, ::WaterBalance, state)

Computes the volumetric water content from current saturation or pressure state.
"""
function watercontent!(sub::SubSurface, water::WaterBalance, state)
    @inbounds for i in eachindex(state.sat)
        state.θsat[i] = maxwater(sub, water, state, i)
        state.θwi[i] = state.sat[i]*state.θsat[i]
        # initially set unfrozen water = water+ice;
        # when heat conduction is included, it should update this based on temperature.
        state.θw[i] = state.θwi[i]
    end
end

"""
    hydraulicconductivity!(sub::SubSurface, water::WaterBalance, state)

Computes hydraulic conductivities for the given subsurface layer and water balance scheme.
"""
function hydraulicconductivity!(sub::SubSurface, water::WaterBalance, state)
    @inbounds for i in eachindex(state.kwc)
        let θsat = Hydrology.maxwater(sub, water, state, i),
            θw = state.θw[i],
            θwi = state.θwi[i],
            kw_sat = kwsat(sub, state, i);
            state.kwc[i] = hydraulicconductivity(water, kw_sat, θw, θwi, θsat)
            if i > 1
                state.kw[i] = min(state.kwc[i-1], state.kwc[i])
            end
        end
    end
    # set hydraulic conductivity at boundaries equal to cell conductivities
    state.kw[1] = state.kwc[1]
    state.kw[end] = state.kwc[end]
end

"""
    limit_upper_flux(water::WaterBalance, jw, θw, θwi, θsat, sat, Δz)

Flux limiter for fluxes coming from "above"; limits `jw` based on the
available water in the current cell as well as free pore space.
"""
function limit_upper_flux(water::WaterBalance, jw, θw, θwi, θsat, sat, Δz)
    # case (i): jw < 0 -> outflow -> limit based on available water
    # case (ii): jw > 0 -> inflow -> limit based on available space
    # case (iii): jw = 0 -> no flow, limit has no effect
    jw = min(max(jw, -θw*Δz), (θsat - θwi)*Δz)
    # influx reduction factors
    r = reduction_factor_in(water, sat)*(jw > zero(jw)) + reduction_factor_out(water, sat)*(jw <= zero(jw))
    return r*jw
end

"""
    limit_lower_flux(water::WaterBalance, jw, θw, θwi, θsat, sat, Δz)

Flux limiter for fluxes coming from "below"; limits `jw` based on the
available water in the current cell as well as free pore space.
"""
function limit_lower_flux(water::WaterBalance, jw, θw, θwi, θsat, sat, Δz)
    # case (i): jw < 0 -> inflow -> limit based on available space
    # case (ii): jw > 0 -> outflow -> limit based on available water
    # case (iii): jw = 0 -> no flow, limit has no effect
    jw = min(max(jw, (θwi - θsat)*Δz), θw*Δz)
    # influx reduction factors
    r = reduction_factor_in(water, sat)*(jw < zero(jw)) + reduction_factor_out(water, sat)*(jw >= zero(jw))
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
    return max(jw, jw_lo, jw_up)*(jw <= zero(jw)) + min(jw, jw_lo, jw_up)*(jw > zero(jw))
end

"""
    balancefluxes!(sub::SubSurface, water::WaterBalance, state)

Sums vertical and evapotranspirative flux components `jw_v` and `jw_ET` and applies `balanceflux` to all grid cell faces.
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
            jw_ET = state.jw_ET[i],
            jw_v = state.jw_v[i],
            jw = (jw_ET + jw_v)*dt;
            state.jw[i] = balanceflux(water, jw, θw_up, θw_lo, θwi_up, θwi_lo, θsat_up, θsat_lo, sat_up, sat_lo, Δz_up, Δz_lo)/dt
        end
    end
    # apply flux limits to uppermost and lowermost edges;
    # this may in some cases be redundant with interact! but is necessary due to the possible addition of ET fluxes
    @inbounds state.jw[1] = limit_upper_flux(water, state.jw[1]*dt, state.θw[1], state.θwi[1], state.θsat[1], state.sat[1], state.Δz[1])/dt
    @inbounds state.jw[end] = limit_lower_flux(water, state.jw[end]*dt, state.θw[end], state.θwi[end], state.θsat[end], state.sat[end], state.Δz[end])/dt
    return nothing
end

"""
    advectiveflux(θw_up, θwi_lo, θsat_lo, θmin, kw)

Computes the advective downward water flux given liquid water content in the upper grid cell (θw_up),
the minimum water content (θmin), and hydraulic conductivity at the boundary (kw).
"""
function advectiveflux(θw_up, θmin, kw)
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
            state.jw_v[i] += advectiveflux(θwᵢ₋₁, θfc, kw)
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

Computes the prognostic time derivative for the water balance, usually based on `dθwi`.
Implementation depends on which water flow scheme is being used.
"""
function waterprognostic!(::SubSurface, ::WaterBalance{<:BucketScheme}, state)
    @inbounds @. state.dsat = state.dθwi / state.θsat
end

# Helper methods
function reduction_factor_out(water::WaterBalance, sat)
    smoothness = water.prop.rf_smoothness
    rf = 1 - exp(-sat/smoothness)
    return rf
end
function reduction_factor_in(water::WaterBalance, sat)
    smoothness = water.prop.rf_smoothness
    rf = 1 - exp((sat-1)/smoothness)
    return rf
end

# reduction_factor(water::WaterBalance, x) = 1 - exp(-water.prop.r_β*(1-x)^2)

"""
    resetfluxes!(::SubSurface, water::WaterBalance, state)

Resets flux terms (`jw` and `dθwi`) for `WaterBalance`.
"""
function CryoGrid.resetfluxes!(::SubSurface, water::WaterBalance, state)
    state.jw .= zero(eltype(state.jw))
    state.jw_v .= zero(eltype(state.jw_v))
    state.jw_ET .= zero(eltype(state.jw_ET))
    state.dθwi .= zero(eltype(state.dθwi))
end

# CryoGrid methods
CryoGrid.variables(water::WaterBalance) = (
    CryoGrid.variables(water.flow)...,
    CryoGrid.variables(water.et)...,
    Diagnostic(:jw, OnGrid(Edges), u"m/s", desc="Total water flux over grid edges."), # water fluxes over grid cell boundaries
    Diagnostic(:jw_v, OnGrid(Edges), u"m/s", desc="Advective/diffusive water flux over grid edges."), # vertical water fluxes over grid cell boundaries
    Diagnostic(:jw_ET, OnGrid(Edges), u"m/s", desc="Water fluxes due to evapotranspiration."), # water fluxes due to evapotranspiration
    Diagnostic(:θwi, OnGrid(Cells), domain=0..1), # total volumetric water+ice content
    Diagnostic(:θw, OnGrid(Cells), domain=0..1), # unfrozen/liquid volumetric water content
    Diagnostic(:θsat, OnGrid(Cells), domain=0..1), # maximum volumetric water content (saturation point)
    Diagnostic(:dθwi, OnGrid(Cells), u"1/s"), # divergence of total water content
    Diagnostic(:kw, OnGrid(Edges), u"m/s", domain=0..Inf), # hydraulic conductivity (edges)
    Diagnostic(:kwc, OnGrid(Cells), u"m/s", domain=0..Inf), # hydraulic conductivity (cells)
)
CryoGrid.variables(::BucketScheme) = (
    Prognostic(:sat, OnGrid(Cells), domain=0..1), # autmoatically generates dsat
)

function CryoGrid.initialcondition!(sub::SubSurface, water::WaterBalance, state)
    CryoGrid.computediagnostic!(sub, water, state)
end

function CryoGrid.computediagnostic!(sub::SubSurface, water::WaterBalance, state)
    watercontent!(sub, water, state)
    hydraulicconductivity!(sub, water, state)
    evapotranspiration!(sub, water, state)
end

function CryoGrid.computefluxes!(sub::SubSurface, water::WaterBalance, state)
    wateradvection!(sub, water, state)
    waterdiffusion!(sub, water, state)
    balancefluxes!(sub, water, state)
    divergence!(state.dθwi, state.jw, Δ(state.grid))
    waterprognostic!(sub, water, state)
end

function CryoGrid.timestep(
    sub::SubSurface,
    water::WaterBalance{TFlow,TET,<:CryoGrid.MaxDelta},
    state
) where {TFlow,TET}
    dtmax = Inf
    @inbounds for i in 1:length(state.sat)
        dt = water.dtlim(state.dθwi[i], state.θwi[i], state.t, zero(state.t), state.θsat[i])
        dt = isfinite(dt) && dt > zero(dt) ? dt : Inf # make sure it's +Inf
        dtmax = min(dtmax, dt)
    end
    return dtmax
end

# No flow case
hydraulicconductivity!(::SubSurface, ::WaterBalance{NoFlow}, state) = nothing

CryoGrid.variables(::NoFlow) = (
    Diagnostic(:sat, OnGrid(Cells), domain=0..1),
    Diagnostic(:θsat, OnGrid(Cells), domain=0..1),
    Diagnostic(:θwi, OnGrid(Cells), domain=0..1),
)

function CryoGrid.initialcondition!(sub::SubSurface, water::WaterBalance{NoFlow}, state)
    @inbounds for i in eachindex(state.sat)
        state.θsat[i] = maxwater(sub, water, state, i)
        state.θwi[i] = state.sat[i]*state.θsat[i]
    end
end

CryoGrid.computefluxes!(::SubSurface, ::WaterBalance{NoFlow}, state) = nothing

CryoGrid.resetfluxes!(::SubSurface, ::WaterBalance{NoFlow}, state) = nothing

CryoGrid.timestep(::SubSurface, ::WaterBalance{NoFlow}, state) = Inf
