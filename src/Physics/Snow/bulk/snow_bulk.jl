threshold(snow::BulkSnowpack) = snow.para.thresh

function snowdensity!(
    snow::BulkSnowpack{<:ConstantDensity},
    mass::SnowMassBalance,
    state
)
    ρsn = snow.para.density.ρsn
    ρw = waterdensity(snow)
    state.ρsn .= ρsn
    state.por .= 1 - ρsn / ρw
    return nothing
end

function Hydrology.watercontent!(snow::BulkSnowpack, ::WaterBalance, state)
    ρw = waterdensity(snow)
    ρsn = snowdensity(snow, state)
    # total water content = snow water + pore water
    @. state.θwi = ρsn / ρw + state.por*state.sat
    # θsat = porespace
    @. state.θsat = state.por
    return nothing
end

# specify single cell (i.e. "bulk") grid
CryoGrid.makegrid(::BulkSnowpack, strategy, bounds) = Grid([bounds[1], bounds[2]])

# Initialization
function CryoGrid.initialcondition!(snow::BulkSnowpack, ::SnowMassBalance, state)
    @. state.Δz = state.dsn = zero(eltype(state.dsn))
    state.sat .= zero(eltype(state.sat))
    return nothing
end

# Events
CryoGrid.events(::BulkSnowpack, ::SnowMassBalance) = (
    ContinuousEvent(:snow_min),
)
# critierion for minimum snow threshold event
function CryoGrid.criterion(
    ::ContinuousEvent{:snow_min},
    snow::BulkSnowpack,
    state,
)
    Δz = getscalar(state.Δz)
    # use threshold adjusted depth as residual;
    # i.e. the event will fire when snow depth crosses this threshold.
    return Δz - threshold(snow)
end

# triggers for minimum snow threshold event
function CryoGrid.trigger!(
    ::ContinuousEvent{:snow_min},
    ::Decreasing,
    snow::BulkSnowpack,
    state
)
    # Case 1: Decreasing snow depth; set everything to zero to remove snowpack
    state.H .= 0.0
    state.T .= 0.0
    state.sat .= 0.0
    state.θwi .= 0.0
    state.swe .= 0.0
    state.dsn .= 0.0
    state.Δz .= 0.0
    return nothing
end
function CryoGrid.trigger!(
    ::ContinuousEvent{:snow_min},
    ::Increasing,
    snow::BulkSnowpack,
    state
)
    # Case 2: Increasing snow depth; initialize temperature and enthalpy state
    # using current upper boundary temperature.
    heatcapacity!(snow, snow.heat, state)
    state.T .= state.T_ub
    state.H .= state.T.*state.C
    state.por .= 1 - getscalar(state.ρsn) / waterdensity(snow)
    state.sat .= zero(eltype(state.sat))
    return nothing
end

# diagnostics

function CryoGrid.computediagnostic!(snow::BulkSnowpack, heat::HeatBalance, state)
    Heat.freezethaw!(snow, heat, state)
    # compute thermal conductivity
    Heat.thermalconductivity!(snow, heat, state)
    if CryoGrid.isactive(snow, state)
        state.k .= state.kc[1]
    end
end

# computefluxes! for free water, enthalpy based HeatBalance on bulk snow layer
function CryoGrid.computefluxes!(
    snow::BulkSnowpack,
    ps::CoupledSnowWaterHeat{TM,TW,<:HeatBalance{FreeWater,<:EnthalpyBased}},
    state
) where {TM,TW}
    mass, water, heat = ps
    computefluxes!(snow, mass, state)
    dsn = getscalar(state.dsn)
    if dsn < snow.para.thresh
        # set divergence to zero if there is no snow
        @. state.dH = zero(eltype(state.H))
    else
        # otherwise call computefluxes! for other processes
        computefluxes!(snow, Coupled(water, heat), state)
    end
    return nothing
end

# ==== Timestep control ==== #

# Snow mass
function CryoGrid.timestep(snow::BulkSnowpack, mass::SnowMassBalance, state)
    dΔz = getscalar(state.dΔz)
    dsn = getscalar(state.dsn)
    thresh = snow.para.thresh
    dtmax = Inf
    if dsn > thresh && dΔz < zero(dΔz)
        dtmax = (dsn - thresh) / abs(dΔz)
    elseif dsn < thresh && dΔz > zero(dΔz)
        dtmax = (thresh - dsn) / dΔz
    end
    return dtmax
end

include("snow_bulk_dynamic.jl")
include("snow_bulk_prescribed.jl")
