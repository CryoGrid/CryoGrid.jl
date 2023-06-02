using CryoGrid.Snow: DynamicSnowMassBalance, LinearAccumulation

"""
    SurfaceWaterBalance{TR,TS} <: BoundaryProcess{Union{WaterBalance, SnowMassBalance}}

The `SurfaceWaterBalance` represents the closure of the water balance at the surface and acts as both
a Neumann-type upper boundary condition for snow and water fluxes as well as an accountant for the
overall water budget.
"""
Base.@kwdef struct SurfaceWaterBalance{TR,TS} <: BoundaryProcess{Union{WaterBalance, SnowMassBalance}}
    rainfall::TR = nothing
    snowfall::TS = nothing
end

"""
Type alias for `WaterHeatBC{TSEB,TSWB} where {TSEB<:SurfaceEnergyBalance,TSWB<:SurfaceWaterBalance}`.
"""
const SurfaceEnergyWaterBalance{TSEB,TSWB} = WaterHeatBC{TSWB,TSEB} where {TSEB<:SurfaceEnergyBalance,TSWB<:SurfaceWaterBalance}

rainfall(swb::SurfaceWaterBalance, t) = 0.0
rainfall(swb::SurfaceWaterBalance{<:Forcing}, t) = swb.rainfall(t)

snowfall(swb::SurfaceWaterBalance, t) = 0.0
snowfall(swb::SurfaceWaterBalance{TR,<:Forcing}, t) where {TR} = swb.snowfall(t)

function infiltrate!(::Top, swb::SurfaceWaterBalance, ::SubSurface, ::WaterBalance, stop, ssub)
    @setscalar stop.jw_infil = min(stop.jw_rain[1], ssub.kw[1])
    ssub.jw[1] += getscalar(stop.jw_infil)
    return nothing
end

function accumulatesnow!(
    ::Top,
    swb::SurfaceWaterBalance,
    ::Snowpack,
    snowmass::DynamicSnowMassBalance{<:LinearAccumulation},
    stop,
    ssnow,
)
    rate_scale = snowmass.accumulation.rate_scale
    ssnow.∂swe∂t[1] += rate_scale*stop.jw_snow[1]
end

function runoff!(::Top, ::SurfaceWaterBalance, state)
    jw_rain = getscalar(state.jw_rain)
    jw_infil = getscalar(state.jw_infil)
    @setscalar state.∂R∂t = (jw_rain - jw_infil)*area(state.grid)
end

function ETflux!(::Top, ::Coupled2{<:SurfaceEnergyBalance,<:SurfaceWaterBalance}, state)
    jw_ET = getscalar(stop.jw_ET)
    @setscalar stop.∂ET∂t = jw_ET*area(stop.grid)
end

CryoGrid.BCKind(::Type{<:SurfaceWaterBalance}) = CryoGrid.Neumann()

CryoGrid.variables(::Top, ::SurfaceWaterBalance) = (
    Prognostic(:R, Scalar, u"m^3", domain=0..Inf),
    Diagnostic(:jw_rain, Scalar, u"m/s", domain=0..Inf),
    Diagnostic(:jw_snow, Scalar, u"m/s", domain=0..Inf),
    Diagnostic(:jw_infil, Scalar, u"m/s", domain=0..Inf),
    Diagnostic(:jw_ET, Scalar, u"m/s", domain=0..Inf),
)

CryoGrid.variables(top::Top, bc::SurfaceEnergyWaterBalance) = (
    CryoGrid.variables(top, bc.water)...,
    CryoGrid.variables(top, bc.heat)...,
    Prognostic(:ET, Scalar, u"m^3"),
)

function CryoGrid.updatestate!(::Top, swb::SurfaceWaterBalance, stop)
    @setscalar stop.jw_snow = snowfall(swb, stop.t)
    @setscalar stop.jw_rain = rainfall(swb, stop.t)
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

function CryoGrid.interact!(
    top::Top,
    swb::SurfaceWaterBalance,
    sub::SubSurface,
    snowmass::SnowMassBalance,
    stop,
    ssub
)
    accumulatesnow!(top, swb, sub, snowmass, stop, ssub)
    return nothing
end
