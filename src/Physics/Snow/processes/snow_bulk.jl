"""
    BulkSnowpack = Snowpack{<:Bulk}

Type alias for Snowpack with `Bulk` parameterization.
"""
const BulkSnowpack = Snowpack{<:Bulk}
# Local alias for HeatConduction Enthalpy type
const Enthalpy = HeatConduction.Enthalpy

threshold(snow::BulkSnowpack) = snow.para.thresh

CryoGrid.callbacks(::BulkSnowpack, ::SnowMassBalance) = (
    SnowFree(),
)
function CryoGrid.criterion(::SnowFree, snow::BulkSnowpack, ::SnowMassBalance, state)
    dsn = getscalar(state.dsn)
    # snow depth shifted down by 1 cm (TODO: should be configurable?);
    # callback will fire when snow depth drops below or rises above this threshold.
    return dsn - one(dsn)*threshold(snow)
end
function CryoGrid.affect!(::SnowFree, snow::BulkSnowpack, ::SnowMassBalance, state)
    @setscalar state.swe = 0.0
    @setscalar state.dsn = 0.0
    @setscalar state.θwi = 0.0
    @setscalar state.H = 0.0
end

# CryoGrid methods
CryoGrid.thickness(::BulkSnowpack, state) = getscalar(state.dsn)
CryoGrid.midpoints(::BulkSnowpack, state) = -getscalar(state.dsn) / 2
# Snow mass balance
CryoGrid.variables(::BulkSnowpack, ::SnowMassBalance{<:Dynamic}) = (
    Prognostic(:swe, Scalar, u"m"),
    Diagnostic(:dsn, Scalar, u"m"),
    Diagnostic(:ρsn, Scalar, u"kg/m^3"),
    Diagnostic(:T_ub, Scalar, u"°C"),
    Diagnostic(:θwi, OnGrid(Cells)),
    Diagnostic(:dθdH, OnGrid(Cells)),
)
CryoGrid.variables(::BulkSnowpack, ::SnowMassBalance{<:Prescribed}) = (
    Diagnostic(:swe, Scalar, u"m"),
    Diagnostic(:dsn, Scalar, u"m"),
    Diagnostic(:ρsn, Scalar, u"kg/m^3"),
    Diagnostic(:T_ub, Scalar, u"°C"),
    Diagnostic(:θwi, OnGrid(Cells)),
)
# diagnostic step for prescribed snow
function CryoGrid.diagnosticstep!(snow::BulkSnowpack, smb::SnowMassBalance{<:Prescribed}, ssnow)
    ρw = 1000.0 # [kg/m^3] TODO: only available in HydroThermalProperties struct... should be moved elsewhere
    @setscalar ssnow.swe = snow_water_equivalent(snow, smb, ssnow)
    @setscalar ssnow.ρsn = snowdensity(snow, smb, ssnow)
    @setscalar ssnow.dsn = getscalar(ssnow.swe)*ρw/getscalar(ssnow.ρsn)
end
# snowfall swe upper boundary
function CryoGrid.interact!(top::Top, bc::Snowfall, snow::BulkSnowpack, smb::SnowMassBalance{<:Dynamic{LinearAccumulation}}, stop, ssnow)
    ssnow.dswe .+= boundaryvalue(bc, top, smb, snow, stop, ssnow)
end
function CryoGrid.prognosticstep!(::BulkSnowpack, smb::SnowMassBalance{<:Dynamic{LinearAccumulation,InstantaneousRunoff}}, ssnow)
    @. state.dswe += state.θwi*state.dθdH*state.dH
end
# Heat dynamics
function CryoGrid.diagnosticstep!(snow::BulkSnowpack, heat::Heat{FreeWater,Enthalpy}, state)
    dsn = getscalar(state.dsn)
    ρsn = getscalar(state.ρsn)
    ρw = heat.prop.ρw
    if dsn > threshold(snow)
        @. state.θwi = ρsn / ρw
        HeatConduction.freezethaw!(snow, heat, state)
    else
        @. state.θwi = 0.0
        @. state.θw = 0.0
        @. state.T = state.T_ub
        @. state.C = heat.prop.ca
    end
    HeatConduction.thermalconductivity!(snow, heat, state)
    @. state.k = state.kc
    return nothing
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
# Special implementation of interact! only for Dirichlet boundary temperature;
# We do this only so that T_ub can be set to the value of the boundary condition
CryoGrid.interact!(top::Top, bc::HeatBC, snow::BulkSnowpack, heat::Heat, stop, ssnow) = CryoGrid.interact!(BoundaryStyle(bc), top, bc, snow, heat, stop, ssnow)
function CryoGrid.interact!(::Dirichlet, top::Top, bc::HeatBC, sub::BulkSnowpack, heat::Heat, stop, ssnow)
    Δk = CryoGrid.thickness(sub, ssnow) # using `thickness` allows for generic layer implementations
    @setscalar ssnow.T_ub = CryoGrid.boundaryvalue(bc, top, heat, sub, stop, ssnow)
    # boundary flux
    @setscalar ssnow.dH_upper = CryoGrid.boundaryflux(bc, top, heat, sub, stop, ssnow)
    @inbounds ssnow.dH[1] += getscalar(ssnow.dH_upper) / Δk[1]
    return nothing # ensure no allocation
end
