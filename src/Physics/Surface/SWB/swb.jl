using CryoGrid.Snow: SnowMassBalance, LinearAccumulation

"""
    SurfaceWaterBalance{TR,TS} <: BoundaryProcess{Union{WaterBalance, SnowMassBalance}}

The `SurfaceWaterBalance` represents the closure of the water balance at the surface and acts as both
a Neumann-type upper boundary condition for snow and water fluxes as well as an accountant for the
water mass balance at the surface.
"""
Base.@kwdef struct SurfaceWaterBalance{TR,TS} <: BoundaryProcess{Union{WaterBalance, SnowMassBalance}}
    rainfall::TR = ConstantForcing(0.0u"m/s")
    snowfall::TS = ConstantForcing(0.0u"m/s")
end

"""
    SurfaceWaterBalance(precip::Forcing{u"m/s",T1}, Tair::Forcing{u"°C",T2}) where {T1,T2}

SWB constructor which derives rainfall and snowfall forcings from total precipitation and air temperature.
"""
function SurfaceWaterBalance(precip::Forcing{u"m/s",T1}, Tair::Forcing{u"°C",T2}) where {T1,T2}
    return SurfaceWaterBalance(
        TimeVaryingForcing(u"m/s", T1, t -> Tair(t) > 0.0 ? precip(t) : 0.0, :rainfall),
        TimeVaryingForcing(u"m/s", T1, t -> Tair(t) <= 0.0 ? precip(t) : 0.0, :snowfall),
    )
end

"""
    SurfaceWaterBalance(forcings::Forcings)

Automatically attempts to extract precipitation forcings from `forcings`.
"""
function SurfaceWaterBalance(forcings::Forcings)
    if hasproperty(forcings, :rainfall) && hasproperty(forcings, :snowfall)
        return SurfaceWaterBalance(forcings.rainfall, forcings.snowfall)
    elseif hasproperty(forcings, :precip) && hasproperty(forcings, :Tair)
        return SurfaceWaterBalance(forcings.precip, forcings.Tair)
    else
        error("forcings must provide either 'rainfall' and 'snowfall' or 'precip' and 'Tair'")
    end
end

function infiltrate!(top::Top, swb::SurfaceWaterBalance, sub::SubSurface, water::WaterBalance, stop, ssub)
    max_infil = Hydrology.maxinfiltration(sub, water, ssub)
    jw_in = min(stop.jw_rain[1], max_infil)
    ssub.jw_v[1] += jw_in
    Hydrology.balancefluxes!(top, swb, sub, water, stop, ssub)
    # set infiltration flux after balancing
    @setscalar stop.jw_infil = ssub.jw[1] - ssub.jw_ET[1]
end

function runoff!(::Top, ::SurfaceWaterBalance, state)
    jw_rain = getscalar(state.jw_rain)
    jw_infil = getscalar(state.jw_infil)
    @setscalar state.drunoff = max(zero(jw_rain), jw_rain - jw_infil)*area(state.grid)
    @setscalar state.dinfil = max(zero(jw_rain), jw_infil)*area(state.grid)
end

CryoGrid.BCKind(::Type{<:SurfaceWaterBalance}) = CryoGrid.Neumann()

CryoGrid.variables(::Top, ::SurfaceWaterBalance) = (
    Prognostic(:runoff, Scalar, u"m^3", domain=0..Inf),
    Prognostic(:infil, Scalar, u"m^3", domain=0..Inf),
    Diagnostic(:jw_rain, Scalar, u"m/s", domain=0..Inf),
    Diagnostic(:jw_snow, Scalar, u"m/s", domain=0..Inf),
    Diagnostic(:jw_infil, Scalar, u"m/s", domain=0..Inf),
    Diagnostic(:jw_ET, Scalar, u"m/s", domain=0..Inf),
)

function CryoGrid.computediagnostic!(::Top, swb::SurfaceWaterBalance, stop)
    @setscalar stop.jw_snow = swb.snowfall(stop.t)
    @setscalar stop.jw_rain = swb.rainfall(stop.t)
end

function CryoGrid.computefluxes!(top::Top, swb::SurfaceWaterBalance, state)
    runoff!(top, swb, state)
end

# interactions

function CryoGrid.interact!(
    top::Top,
    swb::SurfaceWaterBalance,
    sub::SubSurface,
    water::WaterBalance,
    stop,
    ssub
)
    # Case 1: No ET scheme
    ground_ET(::WaterBalance{<:WaterFlow}, ssub) = zero(eltype(ssub.jw))
    # Case 2: with ET
    ground_ET(::WaterBalance{<:WaterFlow,<:Evapotranspiration}, ssub) = ssub.jw_ET[1]
    # flip the sign from the ground ET flux which is positive downward
    @setscalar stop.jw_ET = -ground_ET(water, ssub)
    infiltrate!(top, swb, sub, water, stop, ssub)
    return nothing
end

# surface water + energy balance
include("coupled_sweb.jl")
