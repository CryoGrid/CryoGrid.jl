"""
    Bulk{Tthresh} <: SnowpackParameterization

Simple, bulk ("single layer") snow scheme where snowpack is represented as a single grid cell with homogenous state.
"""
Base.@kwdef struct Bulk{Tthresh} <: SnowpackParameterization
    thresh::Tthresh = 0.02u"m" # snow threshold
end

"""
    BulkSnowpack = Snowpack{<:Bulk}

Type alias for Snowpack with `Bulk` parameterization.
"""
const BulkSnowpack = Snowpack{<:Bulk}
# Local alias for Heat Enthalpy type
const Enthalpy = Heat.Enthalpy

threshold(snow::BulkSnowpack) = snow.para.thresh

CryoGrid.thickness(::BulkSnowpack, state, i::Integer=1) = getscalar(state.dsn)
CryoGrid.midpoint(::BulkSnowpack, state, i::Integer=1) = -getscalar(state.dsn) / 2

# Events
CryoGrid.events(::BulkSnowpack, ::Coupled2{<:SnowMassBalance,<:HeatBalance}) = (
    ContinuousEvent(:snow_min),
)
# critierion for minimum snow threshold event
function CryoGrid.criterion(
    ::ContinuousEvent{:snow_min},
    snow::BulkSnowpack,
    smb::SnowMassBalance,
    state,
)
    # get current snow water equivalent (from forcing or state variable)
    new_swe = getscalar(swe(snow, smb, state))
    # get current snow density
    new_ρsn = getscalar(snowdensity(snow, smb, state))
    # compute actual snow depth/height
    new_dsn = new_swe*snow.prop.ρw/new_ρsn
    # use threshold adjusted depth as residual;
    # i.e. the event will fire when snow depth crosses this threshold.
    return new_dsn - threshold(snow)
end
# triggers for minimum snow threshold event
function CryoGrid.trigger!(
    ::ContinuousEvent{:snow_min},
    ::Decreasing,
    snow::BulkSnowpack,
    ::Coupled2{<:SnowMassBalance,<:HeatBalance},
    state
)
    # Case 1: Decreasing snow depth; set everything to zero to remove snowpack
    state.H .= 0.0
    state.T .= 0.0
    state.θwi .= 0.0
    state.swe .= 0.0
    state.dsn .= 0.0
    return nothing
end
function CryoGrid.trigger!(
    ::ContinuousEvent{:snow_min},
    ::Increasing,
    snow::BulkSnowpack,
    procs::Coupled2{<:SnowMassBalance,<:HeatBalance},
    state
)
    # Case 2: Increasing snow depth; initialize temperature and enthalpy state
    # using current upper boundary temperature.
    _, heat = procs
    θfracs = volumetricfractions(snow, heat, state, 1)
    state.C .= C = Heat.heatcapacity(snow, heat, θfracs...)
    state.T .= state.T_ub
    state.H .= state.T.*C
    return nothing
end
# heat upper boundary (for all bulk implementations)
function CryoGrid.interact!(top::Top, bc::HeatBC, snow::BulkSnowpack, heat::HeatBalance, stop, ssnow)
    CryoGrid.interact!(CryoGrid.BoundaryStyle(bc), top, bc, snow, heat, stop, ssnow)
    return nothing
end
function CryoGrid.interact!(
    ::CryoGrid.Dirichlet,
    top::Top,
    bc::HeatBC,
    snow::BulkSnowpack,
    heat::HeatBalance,
    stop,
    ssnow
)
    @setscalar ssnow.T_ub = CryoGrid.boundaryvalue(bc, top, heat, snow, stop, ssnow)
    if getscalar(ssnow.dsn) < threshold(snow)
        @setscalar ssnow.T = getscalar(ssnow.T_ub)
    end
    # boundary flux
    ssnow.jH[1] += CryoGrid.boundaryflux(bc, top, heat, snow, stop, ssnow)
    return nothing
end

# ==== Dynamic bulk snow scheme ==== #

CryoGrid.variables(snow::BulkSnowpack, smb::DynamicSnowMassBalance) = (
    Prognostic(:swe, Scalar, u"m", domain=0..Inf),
    Diagnostic(:ρsn, Scalar, u"kg/m^3", domain=0..Inf),
    Diagnostic(:θwi, OnGrid(Cells), domain=0..1), 
    CryoGrid.basevariables(snow, smb)...,
)
function CryoGrid.diagnosticstep!(
    snow::BulkSnowpack,
    procs::Coupled(DynamicSnowMassBalance{TAcc,TAbl,TDen}, HeatBalance{FreeWater,<:Enthalpy}),
    state
) where {TAcc,TAbl<:DegreeDayMelt,TDen<:ConstantDensity}
    smb, heat = procs
    ρsn = snow.prop.ρsn_new
    Heat.resetfluxes!(snow, heat, state)
    @setscalar state.θwi = θwi = ρsn / snow.prop.ρw
    @setscalar state.ρsn = ρsn
    dsn = getscalar(state.swe) / θwi
    # only update snowdepth if swe greater than threshold, otherwise, set to zero.
    @setscalar state.dsn = IfElse.ifelse(getscalar(state.swe) >= threshold(snow)*θwi, dsn, zero(dsn))
    # get heat capacity as a function of liquid water content
    f_hc = partial(heatcapacity, Val{:θw}(), snow, heat, state, 1)
    # set temperature and liquid water content according to free water freeze curve,
    # but capping the liquid fraction according to the 'max_unfrozen' parameter.
    max_unfrozen = ablation(smb).max_unfrozen
    θwi_cap = θwi*max_unfrozen
    T, θw, C = Heat.enthalpyinv(heat.freezecurve, f_hc, getscalar(state.H), θwi_cap, heat.prop.L)
    # do not allow temperature to exceed 0°C
    @. state.T = min(T, zero(T))
    @. state.θw = θw
    @. state.C = C
    # compute thermal conductivity
    Heat.thermalconductivity!(snow, heat, state)
    @. state.k = state.kc
    return nothing
end
# snowfall upper boundary
function CryoGrid.interact!(
    top::Top,
    bc::Snowfall,
    snow::BulkSnowpack,
    smb::DynamicSnowMassBalance{<:LinearAccumulation},
    stop,
    ssnow
)
    # upper boundary condition for snow mass balance;
    # apply snowfall to swe
    rate_scale = accumulation(smb).rate_scale
    snowfall_rate = boundaryvalue(bc, top, smb, snow, stop, ssnow)
    @. ssnow.∂swe∂t += rate_scale*snowfall_rate
end
function CryoGrid.prognosticstep!(
    snow::BulkSnowpack,
    smb::DynamicSnowMassBalance{TAcc,TAbl},
    state
) where {TAcc,TAbl<:DegreeDayMelt}
    if getscalar(state.dsn) < threshold(snow)
        # set energy flux to zero if there is no snow
        @. state.∂H∂t = zero(eltype(state.H))
    end
    if getscalar(state.swe) > 0.0 && getscalar(state.T_ub) > 0.0
        ddf = ablation(smb).factor # [m/K/s]
        jH_upper = state.jH[1] # [J/m^3]
        T_ub = getscalar(state.T_ub) # upper boundary temperature
        Tref = 0.0*unit(T_ub) # just in case T_ub has units
        # calculate the melt rate per second via the degree day model
        dmelt = max(ddf*(T_ub-Tref), zero(eltype(state.∂swe∂t))) # [m/s]
        @. state.∂swe∂t += -dmelt
        # set upper heat flux to zero if dmelt > 0;
        # this is due to the energy being (theoretically) "consumed" to melt the snow
        state.jH[1] *= 1 - (dmelt > zero(dmelt))
    end
    return nothing
end

# ==== Prescribed/forced snow scheme ==== #

# Snow mass balance
CryoGrid.variables(snow::BulkSnowpack, smb::PrescribedSnowMassBalance) = (
    Diagnostic(:swe, Scalar, u"m", domain=0..Inf),
    Diagnostic(:ρsn, Scalar, u"kg/m^3", domain=0..Inf),
    Diagnostic(:θwi, OnGrid(Cells), u"kg/m^3", domain=0..1),
    CryoGrid.basevariables(snow, smb)...,
)
CryoGrid.events(::BulkSnowpack, ::Coupled2{<:PrescribedSnowMassBalance,<:HeatBalance}) = (
    ContinuousEvent(:snow_min),
)
function CryoGrid.criterion(
    ::ContinuousEvent{:snow_min},
    snow::BulkSnowpack,
    smb::PrescribedSnowMassBalance,
    state,
)
    new_swe = swe(snow, smb, state)
    new_ρsn = snowdensity(snow, smb, state)
    new_dsn = new_swe*snow.prop.ρw/new_ρsn
    return new_dsn - threshold(snow)
end
function CryoGrid.trigger!(
    ::ContinuousEvent{:snow_min},
    ::Decreasing,
    snow::BulkSnowpack,
    ::Coupled2{<:PrescribedSnowMassBalance,<:HeatBalance},
    state
)
    state.H .= 0.0
end
function CryoGrid.trigger!(
    ::ContinuousEvent{:snow_min},
    ::Increasing,
    snow::BulkSnowpack,
    procs::Coupled2{<:PrescribedSnowMassBalance,<:HeatBalance},
    state
)
    _, heat = procs
    θfracs = volumetricfractions(snow, state, 1)
    C = Heat.heatcapacity(snow, θfracs...)
    state.T .= state.T_ub
    state.H .= state.T.*C
end
function CryoGrid.diagnosticstep!(
    snow::BulkSnowpack,
    procs::Coupled2{<:PrescribedSnowMassBalance,<:HeatBalance{FreeWater,<:Enthalpy}},
    state
)
    smb, heat = procs
    ρw = snow.prop.ρw
    Heat.resetfluxes!(snow, heat, state)
    new_swe = swe(snow, smb, state)
    new_ρsn = snowdensity(snow, smb, state)
    new_dsn = new_swe*ρw/new_ρsn
    @unpack hc_a, kh_a = thermalproperties(snow)
    if new_dsn > threshold(snow)
        # if new snow depth is above threshold, set state variables
        @setscalar state.swe = new_swe
        @setscalar state.ρsn = new_ρsn
        @setscalar state.dsn = new_dsn
        @. state.θwi = new_ρsn / ρw
        Heat.freezethaw!(snow, heat, state)
        # cap temperature at 0°C
        @. state.T = min(state.T, zero(eltype(state.T)))
        Heat.thermalconductivity!(snow, state)
        @. state.k = state.kc
    else
        # otherwise, set to zero
        @setscalar state.swe = 0.0
        @setscalar state.ρsn = 0.0
        @setscalar state.dsn = 0.0
        @setscalar state.θwi = 0.0
        @setscalar state.θw = 0.0
        @setscalar state.C = hc_a
        @setscalar state.kc = kh_a
    end
end
# prognosticstep! for free water, enthalpy based HeatBalance on snow layer
function CryoGrid.prognosticstep!(
    snow::BulkSnowpack,
    ps::Coupled(SnowMassBalance,HeatBalance{FreeWater,<:Enthalpy}),
    state
)
    smb, heat = ps
    prognosticstep!(snow, smb, state)
    dsn = getscalar(state.dsn)
    if dsn < snow.para.thresh
        # set divergence to zero if there is no snow
        @. state.∂H∂t = zero(eltype(state.H))
    else
        # otherwise call prognosticstep! for heat
        prognosticstep!(snow, heat, state)
    end
end
``