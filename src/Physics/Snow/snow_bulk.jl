"""
    Bulk{Tthresh,Theat,Twater} <: SnowpackParameterization

Simple, bulk ("single layer") snow scheme where snowpack is represented as a single grid cell with homogenous state.
"""
Base.@kwdef struct Bulk{Tthresh,Theat,Twater} <: SnowpackParameterization
    thresh::Tthresh = 0.02u"m" # snow threshold
    heat::Theat = ThermalProperties()
    water::Twater = HydraulicProperties()
end

"""
    BulkSnowpack = Snowpack{<:Bulk}

Type alias for Snowpack with `Bulk` parameterization.
"""
const BulkSnowpack = Snowpack{<:Bulk}
# Local alias for Heat Enthalpy type
const Enthalpy = Heat.Enthalpy

threshold(snow::BulkSnowpack) = snow.para.thresh

function snowdensity!(
    snow::BulkSnowpack,
    mass::DynamicSnowMassBalance{TAcc,TAbl,TDen},
    state
) where {TAcc,TAbl,TDen<:ConstantDensity}
    ρsn = snowdensity(snow, mass, state)
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

function Hydrology.watercontent!(snow::BulkSnowpack, ::WaterBalance, state)
    ρw = waterdensity(snow)
    ρsn = snowdensity(snow, snow.mass, state)
    # total water content = snow water + pore water
    @. state.θwi = ρsn / ρw + state.por*state.sat
    # θsat = porespace; note that this 
    state.θsat .= state.por
    return nothing
end

# specify single cell (i.e. "bulk") grid
CryoGrid.makegrid(::BulkSnowpack, strategy, bounds) = Grid([bounds[1], bounds[2]])

# Initialization
function CryoGrid.initialcondition!(::BulkSnowpack, ::SnowMassBalance, state)
    @. state.Δz = state.dsn = zero(eltype(state.dsn))
    @. state.sat = zero(eltype(state.sat))
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
    state.sat .= 0.01
    return nothing
end

# computefluxes! for free water, enthalpy based HeatBalance on snow layer
function CryoGrid.computefluxes!(
    snow::BulkSnowpack,
    ps::Coupled(SnowMassBalance, HeatBalance{FreeWater,<:Enthalpy}),
    state
)
    smb, heat = ps
    computefluxes!(snow, smb, state)
    dsn = getscalar(state.dsn)
    if dsn < snow.para.thresh
        # set divergence to zero if there is no snow
        @. state.∂H∂t = zero(eltype(state.H))
    else
        # otherwise call computefluxes! for heat
        computefluxes!(snow, heat, state)
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

CryoGrid.variables(snow::BulkSnowpack, smb::DynamicSnowMassBalance) = (
    Prognostic(:swe, Scalar, u"m", domain=0..Inf),
    Diagnostic(:ρsn, Scalar, u"kg/m^3", domain=0..Inf),
    Diagnostic(:por, OnGrid(Cells), domain=0..1),
    Diagnostic(:θwi, OnGrid(Cells), domain=0..1),
    snowvariables(snow)...,
)
function CryoGrid.updatestate!(
    snow::BulkSnowpack,
    procs::Coupled(
        DynamicSnowMassBalance{TAcc,TAbl,TDen},
        WaterBalance,
        HeatBalance{FreeWater,<:Enthalpy}
    ),
    state
) where {TAcc,TAbl<:DegreeDayMelt,TDen<:ConstantDensity}
    mass, water, heat = procs
    resetfluxes!(snow, heat, state)
    # update snow density
    snowdensity!(snow, mass, state)
    # update snow depth
    snowdepth!(snow, mass, state)
    # update water content
    Hydrology.watercontent!(snow, water, state)
    # evaluate freezing/thawing processes for snow layer
    Heat.freezethaw!(snow, heat, state)
    # compute thermal conductivity
    Heat.thermalconductivity!(snow, heat, state)
    if CryoGrid.isactive(snow, state)
        state.k .= state.kc[1]
    end
    return nothing
end

function CryoGrid.computefluxes!(
    snow::BulkSnowpack,
    mass::DynamicSnowMassBalance{TAcc,TAbl},
    state
) where {TAcc,TAbl<:DegreeDayMelt}
    if getscalar(state.dsn) < threshold(snow)
        # set energy flux to zero if there is no snow
        @. state.∂H∂t = zero(eltype(state.H))
    else
        ddf = mass.para.ablation.factor # [m/K/s]
        T_ub = getscalar(state.T_ub) # upper boundary temperature
        Tref = 0.0*unit(T_ub) # just in case T_ub has units
        # calculate the melt rate per second via the degree day model
        dmelt = max(ddf*max(T_ub-Tref, zero(T_ub)), zero(eltype(state.∂swe∂t))) # [m/s]
        @. state.∂swe∂t += -dmelt
        # set upper heat flux to zero if dmelt > 0;
        # this is due to the energy being (theoretically) "consumed" to melt the snow
        state.jH[1] *= 1 - (dmelt > zero(dmelt))
    end
    ρsn = snowdensity(snow, mass, state)
    ρw = waterdensity(snow)
    # compute time derivative for moving boundary
    @. state.∂Δz∂t += state.∂swe∂t*ρw/ρsn
    return nothing
end

# ==== Prescribed/forced snow scheme ==== #

# Snow mass balance
CryoGrid.variables(snow::BulkSnowpack, smb::PrescribedSnowMassBalance) = (
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
    smb::PrescribedSnowMassBalance,
    state,
)
    ρw = waterdensity(snow)
    new_swe = swe(snow, smb, state)
    new_ρsn = snowdensity(snow, smb, state)
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
    smb, heat = procs
    ρw = waterdensity(snow)
    resetfluxes!(snow, heat, state)
    new_swe = swe(snow, smb, state)
    new_ρsn = snowdensity(snow, smb, state)
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
