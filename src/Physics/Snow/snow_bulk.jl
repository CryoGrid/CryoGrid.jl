"""
    BulkSnowpack{Tdsn} = Snowpack{<:Bulk{Tdsn}} where {Tdsn}

Type alias for Snowpack with `Bulk` parameterization.
"""
const BulkSnowpack{Tdsn} = Snowpack{<:Bulk{Tdsn}} where {Tdsn}
# Local alias for HeatConduction Enthalpy type
const Enthalpy = HeatConduction.Enthalpy

CryoGrid.variables(snow::BulkSnowpack) = (
    Diagnostic(:θw, OnGrid(Cells)),
    Diagnostic(:T_ub, Scalar, u"°C"),
    Diagnostic(:dsn, Scalar, u"m"),
)

CryoGrid.callbacks(snow::BulkSnowpack, heat::Heat) = (
    SnowFree(),
)
function CryoGrid.criterion(::SnowFree, snow::BulkSnowpack, heat::Heat, state)
    dsn = snowdepth(snow, state)
    # snow depth shifted down by 1 cm (TODO: should be configurable?);
    # callback will fire when snow depth drops below or rises above this threshold.
    return dsn - one(dsn)*snow.para.thresh
end
function CryoGrid.affect!(::SnowFree, snow::BulkSnowpack, heat::Heat, state)
    @. state.θw = 0.0
    @. state.H = 0.0
end

snowdepth(snow::BulkSnowpack{<:Number}, state) = snow.para.dsn
snowdepth(snow::BulkSnowpack{<:Forcing}, state) = snow.para.dsn(state.t)

# HeatConduction methods
function HeatConduction.freezethaw!(snow::Snowpack, heat::Heat{FreeWater,Enthalpy}, state)
    @inbounds for i in 1:length(state.H)
        # liquid water content = (total water content) * (liquid fraction)
        state.θl[i] = HeatConduction.liquidwater(snow, heat, state, i)
        # update heat capacity
        state.C[i] = C = snow.prop.csn
        # enthalpy inverse function
        state.T[i] = HeatConduction.enthalpyinv(snow, heat, state, i)
        # set dHdT (a.k.a dHdT)
        state.dHdT[i] = state.T[i] ≈ 0.0 ? 1e8 : 1/C
    end
    return nothing
end
@inline HeatConduction.thermalconductivities(snow::Snowpack, ::Heat) = (snow.prop.ksn,)
@inline HeatConduction.heatcapacities(snow::Snowpack, ::Heat) = (snow.prop.csn,)

# CryoGrid methods
CryoGrid.volumetricfractions(::BulkSnowpack, ::Heat, state, i) = (1.0,)
CryoGrid.thickness(::BulkSnowpack, state) = getscalar(state.dsn)
CryoGrid.midpoints(::BulkSnowpack, state) = -getscalar(state.dsn) / 2
function CryoGrid.diagnosticstep!(snow::BulkSnowpack, state)
    @setscalar state.dsn = snowdepth(snow, state)
end
function CryoGrid.diagnosticstep!(snow::BulkSnowpack, heat::Heat{FreeWater,Enthalpy}, state)
    dsn = getscalar(state.dsn)
    @. state.θw = snow.prop.ρsn / heat.prop.ρw
    if dsn > snow.para.thresh
        @. state.θw = snow.prop.ρsn / heat.prop.ρw
        HeatConduction.freezethaw!(snow, heat, state)
    else
        @. state.θw = 0.0
        @. state.T = state.T_ub
        @. state.C = heat.prop.ca
    end
    HeatConduction.thermalconductivity!(snow, heat, state)
    @. state.k = state.kc
    return nothing
end
function CryoGrid.prognosticstep!(snow::BulkSnowpack, heat::Heat{FreeWater,Enthalpy}, state)
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
