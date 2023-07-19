"""
    Bulk{Tden,Tthresh,Theat,Twater} <: SnowpackParameterization

Simple, bulk ("single layer") snow scheme where snowpack is represented as a single grid cell with homogenous state.
"""
Base.@kwdef struct Bulk{Tden<:SnowDensityScheme,Tthresh,Theat,Twater} <: SnowpackParameterization
    thresh::Tthresh = 0.02u"m" # snow threshold
    density::Tden = ConstantDensity() # snow density
    heat::Theat = ThermalProperties() # thermal properties
    water::Twater = HydraulicProperties(kw_sat=1e-4) # hydraulic properties
end

"""
    BulkSnowpack = Snowpack{<:Bulk}

Type alias for Snowpack with `Bulk` parameterization.
"""
const BulkSnowpack{T} = Snowpack{<:Bulk{T}} where {T}
# Local alias for Heat Enthalpy type
const Enthalpy = Heat.Enthalpy

threshold(snow::BulkSnowpack) = snow.para.thresh

function snowdensity!(
    snow::BulkSnowpack{<:ConstantDensity},
    mass::DynamicSnowMassBalance,
    state
)
    ρsn = snow.para.density.ρsn
    ρw = waterdensity(snow)
    state.ρsn .= ρsn
    state.por .= 1 - ρsn / ρw
    return nothing
end

function snowdepth!(
    ::BulkSnowpack,
    ::DynamicSnowMassBalance,
    state
)
    @setscalar state.dsn = getscalar(state.Δz)
end

# implement ablation! for DegreeDayMelt
function ablation!(
    ::Top,
    ::SnowBC,
    snow::BulkSnowpack,
    mass::DynamicSnowMassBalance{TAcc,<:DegreeDayMelt},
    stop,
    ssnow,
) where {TAcc}
    if isactive(snow, ssnow)
        T_ub = getscalar(ssnow.T_ub) # upper boundary temperature
        dmelt = calculate_degree_day_snow_melt(mass.ablation, T_ub)
        dmelt = min(dmelt, getscalar(ssnow.swe))
        # swe flux
        @. ssnow.∂swe∂t -= dmelt
        # thickness flux
        por = getscalar(ssnow.por)
        θis = 1 - por # solid ice
        Δdsn = -dmelt / θis
        @. ssnow.∂Δz∂t += Δdsn
        # add water flux due to melt
        sat = getscalar(ssnow.sat)
        ssnow.jw[1] += dmelt - Δdsn*por*sat
    end
end

# simple linear accumulation scheme for bulk snow
function accumulate!(
    ::Top,
    snowbc::SnowBC,
    snowpack::BulkSnowpack,
    mass::DynamicSnowMassBalance{<:LinearAccumulation},
    stop,
    ssnow,
)
    rate_scale = mass.accumulation.rate_scale
    jw_snow = snowfall(snowbc, stop)
    Δswe = rate_scale*jw_snow
    @. ssnow.∂swe∂t += Δswe
    por = getscalar(ssnow.por)
    θis = 1 - por # solid ice
    Δdsn = Δswe/ θis
    @. ssnow.∂Δz∂t += Δdsn
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
    heat = snow.heat
    θfracs = volumetricfractions(snow, state, 1)
    state.C .= C = Heat.heatcapacity(snow, heat, θfracs...)
    state.T .= state.T_ub
    state.H .= state.T.*C
    state.por .= 1 - getscalar(state.ρsn) / waterdensity(snow)
    state.sat .= zero(eltype(state.sat))
    return nothing
end

# computefluxes! for free water, enthalpy based HeatBalance on bulk snow layer
function CryoGrid.computefluxes!(
    snow::BulkSnowpack,
    ps::Coupled(SnowMassBalance, HeatBalance{FreeWater,<:Enthalpy}),
    state
)
    mass, heat = ps
    computefluxes!(snow, mass, state)
    dsn = getscalar(state.dsn)
    if dsn < snow.para.thresh
        # set divergence to zero if there is no snow
        @. state.∂H∂t = zero(eltype(state.H))
    else
        computefluxes!(snow, heat, state)
    end
    return nothing
end
function CryoGrid.computefluxes!(
    snow::BulkSnowpack,
    ps::Coupled(SnowMassBalance, WaterBalance, HeatBalance{FreeWater,<:Enthalpy}),
    state
)
    mass, water, heat = ps
    computefluxes!(snow, mass, state)
    dsn = getscalar(state.dsn)
    if dsn < snow.para.thresh
        # set divergence to zero if there is no snow
        @. state.∂H∂t = zero(eltype(state.H))
    else
        # otherwise call computefluxes! for other processes
        computefluxes!(snow, Coupled(water, heat), state)
    end
    return nothing
end

# Timestep control
CryoGrid.timestep(::Snowpack, heat::HeatBalance{<:FreeWater,THeatOp,<:CryoGrid.CFL}, state) where {THeatOp} = error("CFL is not supported on snow layer")
function CryoGrid.timestep(snow::Snowpack, heat::HeatBalance{<:FreeWater,THeatOp,<:CryoGrid.MaxDelta}, state) where {THeatOp}
    Δx = Δ(state.grid)
    dtmax = Inf
    if getscalar(state.dsn) > snow.para.thresh
        @inbounds for i in eachindex(Δx)
            dtmax = min(dtmax, heat.dtlim(state.∂H∂t[i], state.H[i], state.t))
        end
        dtmax = isfinite(dtmax) && dtmax > 0 ? dtmax : Inf
    end
    return dtmax
end

# ==== Dynamic bulk snow scheme ==== #

CryoGrid.variables(snow::BulkSnowpack, ::DynamicSnowMassBalance) = (
    Prognostic(:swe, Scalar, u"m", domain=0..Inf),
    Diagnostic(:ρsn, Scalar, u"kg/m^3", domain=0..Inf),
    Diagnostic(:por, OnGrid(Cells), domain=0..1),
    Diagnostic(:θwi, OnGrid(Cells), domain=0..1),
    snowvariables(snow)...,
)
function CryoGrid.updatestate!(
    snow::BulkSnowpack{<:ConstantDensity},
    procs::Coupled(
        DynamicSnowMassBalance{TAcc,<:DegreeDayMelt},
        WaterBalance,
        HeatBalance{FreeWater,<:Enthalpy}
    ),
    state
) where {TAcc}
    mass, water, heat = procs
    # update snow density
    snowdensity!(snow, mass, state)
    # update snow depth
    snowdepth!(snow, mass, state)
    # update water content
    updatestate!(snow, water, state)
    # evaluate freezing/thawing processes for snow layer
    Heat.freezethaw!(snow, heat, state)
    # compute thermal conductivity
    Heat.thermalconductivity!(snow, heat, state)
    if CryoGrid.isactive(snow, state)
        state.k .= state.kc[1]
    end
    return nothing
end

# ==== Prescribed/forced snow scheme ==== #

# Snow mass balance
CryoGrid.variables(snow::BulkSnowpack, ::PrescribedSnowMassBalance) = (
    Diagnostic(:swe, Scalar, u"m", domain=0..Inf),
    Diagnostic(:ρsn, Scalar, u"kg/m^3", domain=0..Inf),
    Diagnostic(:θwi, OnGrid(Cells), u"kg/m^3", domain=0..1),
    snowvariables(snow)...,
)
CryoGrid.events(::BulkSnowpack, ::Coupled(PrescribedSnowMassBalance, HeatBalance)) = (
    ContinuousEvent(:snow_min),
)
function CryoGrid.criterion(
    ::ContinuousEvent{:snow_min},
    snow::BulkSnowpack,
    mass::PrescribedSnowMassBalance,
    state,
)
    ρw = waterdensity(snow)
    new_swe = swe(snow, mass, state)
    new_ρsn = snowdensity(snow, state)
    new_dsn = new_swe*ρw/new_ρsn
    return new_dsn - threshold(snow)
end
function CryoGrid.trigger!(
    ::ContinuousEvent{:snow_min},
    ::Decreasing,
    snow::BulkSnowpack,
    ::Coupled(PrescribedSnowMassBalance,HeatBalance),
    state
)
    state.H .= 0.0
    return nothing
end
function CryoGrid.trigger!(
    ::ContinuousEvent{:snow_min},
    ::Increasing,
    snow::BulkSnowpack,
    procs::Coupled(PrescribedSnowMassBalance, HeatBalance),
    state
)
    _, heat = procs
    θfracs = volumetricfractions(snow, state, 1)
    C = Heat.heatcapacity(snow, heat, θfracs...)
    state.T .= state.T_ub
    state.H .= state.T.*C
    return nothing
end
function CryoGrid.updatestate!(
    snow::BulkSnowpack,
    procs::Coupled(PrescribedSnowMassBalance,HeatBalance{FreeWater,<:Enthalpy}),
    state
)
    mass, heat = procs
    ρw = waterdensity(snow)
    new_swe = swe(snow, mass, state)
    new_ρsn = snowdensity(snow, state)
    new_dsn = new_swe*ρw/new_ρsn
    @setscalar state.Δz = new_dsn
    @unpack ch_a, kh_a = thermalproperties(snow)
    if new_dsn > threshold(snow)
        # if new snow depth is above threshold, set state variables
        @setscalar state.swe = new_swe
        @setscalar state.ρsn = new_ρsn
        @setscalar state.dsn = new_dsn
        @. state.θwi = new_ρsn / ρw
        Heat.freezethaw!(snow, heat, state)
        # cap temperature at 0°C
        @. state.T = min(state.T, zero(eltype(state.T)))
        Heat.thermalconductivity!(snow, heat, state)
        @. state.k = state.kc
    else
        # otherwise, set to zero
        @setscalar state.swe = 0.0
        @setscalar state.ρsn = 0.0
        @setscalar state.dsn = 0.0
        state.θwi .= 0.0
        state.θw .= 0.0
        state.C .= ch_a
        state.kc .= kh_a
    end
    return nothing
end
