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
# Local alias for HeatConduction Enthalpy type
const Enthalpy = HeatConduction.Enthalpy

threshold(snow::BulkSnowpack) = snow.para.thresh

CryoGrid.thickness(::BulkSnowpack, state, i::Integer=1) = getscalar(state.dsn)
CryoGrid.midpoint(::BulkSnowpack, state, i::Integer=1) = -getscalar(state.dsn) / 2

# Events
CryoGrid.events(::BulkSnowpack, ::Coupled2{<:SnowMassBalance,<:Heat}) = (
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
    ::Coupled2{<:SnowMassBalance,<:Heat},
    state
)
    # Case 1: Decreasing snow depth; set everything to zero to remove snowpack
    state.H .= 0.0
    state.T .= 0.0
    state.θwi .= 0.0
    state.swe .= 0.0
    state.dsn .= 0.0
end
function CryoGrid.trigger!(
    ::ContinuousEvent{:snow_min},
    ::Increasing,
    snow::BulkSnowpack,
    procs::Coupled2{<:SnowMassBalance,<:Heat},
    state
)
    # Case 2: Increasing snow depth; initialize temperature and enthalpy state
    # using current upper boundary temperature.
    _, heat = procs
    θfracs = volumetricfractions(snow, heat, state, 1)
    state.C .= C = HeatConduction.heatcapacity(snow, heat, θfracs...)
    state.T .= state.T_ub
    state.H .= state.T.*C
end
# heat upper boundary (for all bulk implementations)
CryoGrid.interact!(top::Top, bc::HeatBC, snow::BulkSnowpack, heat::Heat, stop, ssnow) = CryoGrid.interact!(CryoGrid.BoundaryStyle(bc), top, bc, snow, heat, stop, ssnow)
function CryoGrid.interact!(::CryoGrid.Dirichlet, top::Top, bc::HeatBC, snow::BulkSnowpack, heat::Heat, stop, ssnow)
    Δk = CryoGrid.thickness(snow, ssnow) # using `thickness` allows for generic layer implementations
    @setscalar ssnow.T_ub = CryoGrid.boundaryvalue(bc, top, heat, snow, stop, ssnow)
    if getscalar(ssnow.dsn) < threshold(snow)
        @setscalar ssnow.T = getscalar(ssnow.T_ub)
    end
    # boundary flux
    @setscalar ssnow.dH_upper = CryoGrid.boundaryflux(bc, top, heat, snow, stop, ssnow)
    @inbounds ssnow.dH[1] += getscalar(ssnow.dH_upper) / Δk[1]
    return nothing # ensure no allocation
end

# ==== Dynamic bulk snow scheme ==== #

CryoGrid.variables(snow::BulkSnowpack, smb::DynamicSnowMassBalance) = (
    Prognostic(:swe, Scalar, u"m", domain=0..Inf),
    Diagnostic(:ρsn, Scalar, u"kg/m^3", domain=0..Inf),
    Diagnostic(:θwi, OnGrid(Cells), u"kg/m^3", domain=0..1), 
    snowvariables(snow, smb)...,
)
function CryoGrid.diagnosticstep!(
    snow::BulkSnowpack,
    procs::Coupled2{<:DynamicSnowMassBalance{TAcc,TAbl,TDen},<:Heat{FreeWater,Enthalpy}},
    state
) where {TAcc,TAbl<:DegreeDayMelt,TDen<:ConstantDensity}
    smb, heat = procs
    ρsn = snow.prop.ρsn_new
    @setscalar state.θwi = θwi = ρsn / snow.prop.ρw
    @setscalar state.ρsn = ρsn
    dsn = getscalar(state.swe) / θwi
    # only update snowdepth if swe greater than threshold, otherwise, set to zero.
    @setscalar state.dsn = IfElse.ifelse(getscalar(state.swe) >= threshold(snow)*θwi, dsn, zero(dsn))
    # get heat capacity as a function of liquid water content
    f_hc = partial(heatcapacity, liquidwater, snow, heat, state, 1)
    # set temperature and liquid water content according to free water freeze curve,
    # but capping the liquid fraction according to the 'max_unfrozen' parameter.
    max_unfrozen = ablation(smb).max_unfrozen
    θwi_cap = θwi*max_unfrozen
    T, θw, C = HeatConduction.enthalpyinv(heat.freezecurve, f_hc, getscalar(state.H), heat.prop.L, θwi_cap)
    # do not allow temperature to exceed 0°C
    @. state.T = min(T, zero(T))
    @. state.θw = θw
    @. state.C = C
    # compute thermal conductivity
    HeatConduction.thermalconductivity!(snow, heat, state)
    @. state.k = state.kc
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
    @. ssnow.dswe += rate_scale*snowfall_rate
end
function CryoGrid.prognosticstep!(
    snow::BulkSnowpack,
    smb::DynamicSnowMassBalance{TAcc,TAbl},
    state
) where {TAcc,TAbl<:DegreeDayMelt}
    if getscalar(state.dsn) < threshold(snow)
        # set energy flux to zero if there is no snow
        @. state.dH = zero(eltype(state.H))
    end
    if getscalar(state.swe) > 0.0 && getscalar(state.T_ub) > 0.0
        ddf = ablation(smb).factor # [m/K/s]
        dH_upper = getscalar(state.dH_upper) # [J/m^3]
        T_ub = getscalar(state.T_ub) # upper boundary temperature
        Tref = 0.0*unit(T_ub) # just in case T_ub has units
        # calculate the melt rate per second via the degree day model
        dmelt = max(ddf*(T_ub-Tref), zero(dH_upper)) # [m/s]
        @. state.dswe += -dmelt
        # remove upper heat flux from dH if dmelt > 0;
        # this is due to the energy being "consumed" to melt the snow
        @. state.dH += -dH_upper*(dmelt > zero(dmelt))
    end
    return nothing
end

# ==== Prescribed/forced snow scheme ==== #

# Snow mass balance
CryoGrid.variables(snow::BulkSnowpack, smb::PrescribedSnowMassBalance) = (
    Diagnostic(:swe, Scalar, u"m", domain=0..Inf),
    Diagnostic(:ρsn, Scalar, u"kg/m^3", domain=0..Inf),
    Diagnostic(:θwi, OnGrid(Cells), u"kg/m^3", domain=0..1),
    snowvariables(snow, smb)...,
)
CryoGrid.events(::BulkSnowpack, ::Coupled2{<:PrescribedSnowMassBalance,<:Heat}) = (
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
    ::Coupled2{<:PrescribedSnowMassBalance,<:Heat},
    state
)
    state.H .= 0.0
end
function CryoGrid.trigger!(
    ::ContinuousEvent{:snow_min},
    ::Increasing,
    snow::BulkSnowpack,
    procs::Coupled2{<:PrescribedSnowMassBalance,<:Heat},
    state
)
    _, heat = procs
    θfracs = volumetricfractions(snow, heat, state, 1)
    C = HeatConduction.heatcapacity(snow, heat, θfracs...)
    state.T .= state.T_ub
    state.H .= state.T.*C
end
function CryoGrid.diagnosticstep!(
    snow::BulkSnowpack,
    procs::Coupled2{<:PrescribedSnowMassBalance,<:Heat{FreeWater,Enthalpy}},
    state
)
    smb, heat = procs
    ρw = snow.prop.ρw
    new_swe = swe(snow, smb, state)
    new_ρsn = snowdensity(snow, smb, state)
    new_dsn = new_swe*ρw/new_ρsn
    ρw = heat.prop.consts.ρw
    if new_dsn > threshold(snow)
        # if new snow depth is above threshold, set state variables
        @setscalar state.swe = new_swe
        @setscalar state.ρsn = new_ρsn
        @setscalar state.dsn = new_dsn
        @. state.θwi = new_ρsn / ρw
        HeatConduction.freezethaw!(snow, heat, state)
        # cap temperature at 0°C
        @. state.T = min(state.T, zero(eltype(state.T)))
        HeatConduction.thermalconductivity!(snow, heat, state)
        @. state.k = state.kc
    else
        # otherwise, set to zero
        @setscalar state.swe = 0.0
        @setscalar state.ρsn = 0.0
        @setscalar state.dsn = 0.0
        @setscalar state.θwi = 0.0
        @setscalar state.θw = 0.0
        @setscalar state.C = heat.prop.ca
        @setscalar state.kc = heat.prop.ka
    end
end
# override prognosticstep! for incompatible heat types to prevent incorrect usage;
# this will actually result in an ambiguous dispatch error before this error is thrown.
CryoGrid.prognosticstep!(snow::BulkSnowpack, heat::Heat, state) = error("prognosticstep! not implemented for $(typeof(heat)) on $(typeof(snow))")
# prognosticstep! for free water, enthalpy based Heat on snow layer
function CryoGrid.prognosticstep!(snow::BulkSnowpack, ::Heat{FreeWater,Enthalpy}, state)
    dsn = getscalar(state.dsn)
    if dsn < snow.para.thresh
        # set energy flux to zero if there is no snow
        @. state.dH = zero(eltype(state.H))
    end
end
