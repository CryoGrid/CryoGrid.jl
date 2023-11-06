# ==== Prescribed/forced snow scheme ==== #

function compute_diagnostic_prescribed_snow_heat!(
    snow::Snowpack{<:Bulk,<:PrescribedSnowMassBalance},
    heat::HeatBalance{FreeWater,<:EnthalpyBased},
    state
)
    mass = snow.mass
    dsn = getscalar(snowdepth(snow, mass, state))
    @unpack ch_a, kh_a = thermalproperties(snow)
    if dsn > threshold(snow)
        Heat.freezethaw!(snow, heat, state)
        # cap temperature at 0°C
        @. state.T = min(state.T, zero(eltype(state.T)))
        Heat.thermalconductivity!(snow, heat, state)
        @. state.k = state.kc
    else
        state.θw .= 0.0
        state.C .= ch_a
        state.kc .= kh_a
    end
end

# Snow methods

function CryoGrid.computediagnostic!(
    snow::BulkSnowpack,
    mass::PrescribedSWE,
    state,
)
    # update snow density
    snowdensity!(snow, mass, state)
    # get snow and water density
    ρsn = getscalar(snowdensity(snow, state))
    ρw = waterdensity(snow)
    # set current swe based on forcing
    state.swe .= swe = snowwater(snow, mass, state)
    # derive snow depth from swe and density
    state.dsn .= swe*ρw / ρsn
end

function CryoGrid.computediagnostic!(
    snow::BulkSnowpack,
    mass::PrescribedSnowDepth,
    state,
)
    # update snow density
    snowdensity!(snow, mass, state)
    # get snow and water density
    ρsn = getscalar(snowdensity(snow, state))
    ρw = waterdensity(snow)
    # set current snow depth based on forcing
    state.dsn .= dsn = snowdepth(snow, mass, state)
    # derive swe from depth and density
    state.swe .= dsn*ρsn / ρw
end

# CryoGrid methods

CryoGrid.variables(snow::BulkSnowpack, ::PrescribedSnowMassBalance) = (
    Diagnostic(:swe, Scalar, u"m", domain=0..Inf),
    Diagnostic(:ρsn, Scalar, u"kg/m^3", domain=0..Inf),
    Diagnostic(:θwi, OnGrid(Cells), u"m^3/m^3", domain=0..1),
    snowvariables(snow)...,
)

CryoGrid.events(::BulkSnowpack, ::Coupled(PrescribedSnowMassBalance, HeatBalance)) = (
    ContinuousEvent(:snow_min),
)

function CryoGrid.criterion(
    ::ContinuousEvent{:snow_min},
    snow::BulkSnowpack,
    state,
)
    dsn = getscalar(snowdepth(snow, snow.mass, state))
    return dsn - threshold(snow)
end

function CryoGrid.trigger!(
    ::ContinuousEvent{:snow_min},
    ::Decreasing,
    snow::BulkSnowpack,
    state
)
    state.H .= 0.0
    return nothing
end

function CryoGrid.trigger!(
    ::ContinuousEvent{:snow_min},
    ::Increasing,
    snow::BulkSnowpack,
    state
)
    c_snow = heatcapacity(snow, snow.heat, state, 1)
    state.T .= state.T_ub
    state.H .= state.T.*c_snow
    return nothing
end

function CryoGrid.computediagnostic!(
    snow::BulkSnowpack,
    procs::CoupledSnowWaterHeat{
        <:PrescribedSnowMassBalance,
        <:WaterBalance,
        <:HeatBalance{FreeWater,<:EnthalpyBased}
    },
    state
)
    mass, water, heat = procs
    computediagnostic!(snow, mass, state)
    computediagnostic!(snow, water, state)
    # special implementation of heat diagnostics that accounts for presence of snow
    compute_diagnostic_prescribed_snow_heat!(snow, heat, state)
    return nothing
end
