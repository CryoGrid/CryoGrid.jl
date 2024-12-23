"""
    Bulk{Tden,Tthresh,Theat,Twater} <: SnowpackParameterization

Simple, bulk ("single layer") snow scheme where snowpack is represented as a single grid cell with homogenous state.
"""
Base.@kwdef struct Bulk{Tden<:SnowDensityScheme,Tthresh,Theat,Twater} <: SnowpackParameterization
    thresh::Tthresh = 0.005u"m" # snow threshold
    density::Tden = ConstantDensity() # snow density
    heat::Theat = SnowThermalProperties() # thermal properties
    water::Twater = HydraulicProperties(kw_sat=1e-4) # hydraulic properties
end

"""
    BulkSnowpack = Snowpack{<:Bulk}

Type alias for Snowpack with `Bulk` parameterization.
"""
const BulkSnowpack{TD} = Snowpack{<:Bulk{TD}} where {TD}

"""
    threshold(snow::BulkSnowpack)

Retrieves the minimum snow threshold for the bulk snow scheme.
"""
threshold(snow::BulkSnowpack) = snow.para.thresh

# do nothing for constant density scheme
compaction!(::Top, ::SnowBC, ::BulkSnowpack{<:ConstantDensity}, ::SnowMassBalance, stop, ssnow) = nothing

function snowdensity!(
    snow::BulkSnowpack{<:ConstantDensity},
    mass::SnowMassBalance,
    state
)
    ρsn = snow.para.density.ρsn
    ρw = waterdensity(snow)
    state.ρsn .= ρsn
    state.θsat .= 1 - ρsn / ρw
    return nothing
end

# implement ablation! for DegreeDayMelt
function ablation!(
    ::Top,
    ::SnowBC,
    snow::BulkSnowpack,
    mass::SnowMassBalance{TAcc,<:DegreeDayMelt},
    stop,
    ssnow,
) where {TAcc}
    swe = getscalar(ssnow.swe)
    T_ub = getscalar(stop.T_ub) # upper boundary temperature
    dd_melt = calculate_degree_day_snow_melt(mass.ablation, T_ub)
    dmelt = min(dd_melt, max(swe, zero(swe)))
    # swe flux
    @. ssnow.dswe -= dmelt
    θsat = getscalar(ssnow.θsat)
    θis = 1 - θsat # solid ice
    Δdsn = -dmelt / θis
    # add water flux due to melt if and only if snowpack is "active"
    sat = getscalar(ssnow.sat)
    water_flux_in = (dmelt + Δdsn*θsat*sat)*isactive(snow, ssnow)
    ssnow.jw[1] += water_flux_in
    return nothing
end

# simple linear accumulation scheme for bulk snow
function accumulation!(
    ::Top,
    snowbc::SnowBC,
    snowpack::BulkSnowpack,
    mass::SnowMassBalance{<:LinearAccumulation},
    stop,
    ssnow,
)
    # get scaling factor(s)
    rate_scale = mass.accumulation.rate_scale
    jw_snow = snowfall(snowbc, stop)[1]
    Δswe = rate_scale*jw_snow
    @. ssnow.dswe += Δswe
    # θsat = getscalar(ssnow.θsat)
    # θis = 1 - θsat # solid ice
    # Δdsn = Δswe / θis
end

function Hydrology.watercontent!(snow::BulkSnowpack, ::WaterBalance, state)
    ρw = waterdensity(snow)
    ρsn = snowdensity(snow, state)
    θis = ρsn / ρw
    # note that if the snow density changes, this will actually
    # water will be lost; ideally the saturation should also change
    state.θsat .= 1 .- ρsn / ρw
    # total water content = snow water + pore water
    state.θwi .= θis + state.θsat.*state.sat
    return nothing
end

CryoGrid.thickness(::BulkSnowpack, state, i::Integer=1) = getscalar(state.dsn)

CryoGrid.midpoint(::BulkSnowpack, state, i::Integer=1) = getscalar(cells(state.grid))

CryoGrid.variables(snow::BulkSnowpack, ::SnowMassBalance) = (
    Prognostic(:swe, Scalar, u"m", domain=0..Inf),
    Diagnostic(:ρsn, Scalar, u"kg/m^3", domain=0..Inf),
    Diagnostic(:θwi, OnGrid(Cells), domain=0..1),
    snowvariables(snow)...,
)

# specify single cell (i.e. "bulk") grid
CryoGrid.makegrid(::BulkSnowpack, strategy, bounds) = Grid([bounds[1], bounds[2]])

# Initialization
function CryoGrid.initialcondition!(snow::BulkSnowpack, ::SnowMassBalance, state)
    state.T .= min(state.T_ub, zero(state.T_ub))
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
    Δz = getscalar(Δ(state.grid))
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
    cellthick = Δ(state.grid)
    # Case 1: Decreasing snow depth; set everything to zero to remove snowpack
    state.H .= 0.0
    state.T .= 0.0
    state.sat .= 0.0
    state.θwi .= 0.0
    state.swe .= 0.0
    state.dsn .= 0.0
    cellthick .= 0.0
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
    state.θsat .= 1 - getscalar(state.ρsn) / waterdensity(snow)
    state.sat .= zero(eltype(state.sat))
    return nothing
end

# mass balance fluxes are handled in interact! for bulk snowpack
CryoGrid.computeprognostic!(snow::BulkSnowpack, mass::SnowMassBalance, state) = nothing

# computeprognostic! for free water, enthalpy based HeatBalance on bulk snow layer
function CryoGrid.computeprognostic!(
    snow::BulkSnowpack,
    ps::CoupledSnowWaterHeat{TM,TW,TH},
    state
) where {TM,TW,TH<:HeatBalance{<:EnthalpyBased}}
    mass, water, heat = ps
    computeprognostic!(snow, mass, state)
    dsn = getscalar(state.dsn)
    if dsn < snow.para.thresh
        # set fluxes to zero if there is no snow
        state.dH .= 0.0
        if !isa(water.flow, NoFlow)
            state.dsat .= 0.0
        end
    else
        # otherwise call computeprognostic! for coupled water/heat
        computeprognostic!(snow, Coupled(water, heat), state)
    end
    return nothing
end

# diagnostic
function CryoGrid.computediagnostic!(snow::BulkSnowpack, heat::HeatBalance, state)
    Heat.freezethaw!(freezecurve(snow), snow, heat, state)
    # compute thermal conductivity
    Heat.thermalconductivity!(snow, state)
    if CryoGrid.isactive(snow, state)
        state.k .= state.kc[1]
    end
end

function CryoGrid.diagnosticstep!(snow::BulkSnowpack, state)
    # force swe to be >= 0
    swe = getscalar(state.swe)
    dsn = getscalar(state.dsn)
    cellthick = Δ(state.grid)
    if swe < zero(swe)
        state.swe .= zero(swe)
        cellthick .= state.dsn .= zero(dsn)
        return true
    end
    return false
end

# ==== Timestep control ==== #

# Snow mass
function CryoGrid.timestep(snow::BulkSnowpack, mass::SnowMassBalance, state)
    ρw = waterdensity(snow)
    ρsn = getscalar(snowdensity(snow, state))
    dΔz = getscalar(state.dswe)*ρsn/ρw
    dsn = getscalar(state.dsn)
    thresh = snow.para.thresh
    dtmax = Inf
    if dsn > thresh && dΔz < zero(dΔz)
        dtmax = (dsn - thresh) / abs(dΔz)
    end
    # take min of threshold dtmax and dt limiter
    # dtmax = mass.dtlim(dΔz, dsn, state.t, 0.0, Inf)
    return dtmax
end
