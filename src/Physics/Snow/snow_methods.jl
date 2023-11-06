### Snow methods ###
# Mandatory
"""
    snowdensity!(::Snowpack, ::SnowMassBalance, state)

Computes the current snow density (if necessary) and stores the result in `state.ρsn`.
"""
snowdensity!(::Snowpack, ::SnowMassBalance, state) = error("not implemented")

"""
    ablation!(::Top, ::SnowBC, ::Snowpack, ::SnowMassBalance, stop, ssnow)

Computes snow mass balance fluxes due to ablation (e.g. snow melt).
"""
ablation!(::Top, ::SnowBC, ::Snowpack, ::DynamicSnowMassBalance, stop, ssnow) = error("not implemented")

"""
    accumulation!(::Top, ::SnowBC, ::Snowpack, ::SnowMassBalance, stop, ssnow) 

Computes snow mass balance fluxes due to accumulation (e.g. snowfall).
"""
accumulation!(::Top, ::SnowBC, ::Snowpack, ::DynamicSnowMassBalance, stop, ssnow) = error("not implemented")

# Optional (w/ default implementations)

"""
    threshold(::Snowpack)

Retrieves the snow cover threshold for this `Snowpack` layer to become active.
"""
threshold(::Snowpack) = 0.01 # meters

"""
    snowwater(::Snowpack, ::SnowMassBalance, state)

Retrieve the current snow water equivalent of the snowpack.
"""
snowwater(::Snowpack, ::SnowMassBalance, state) = state.swe

"""
    snowdensity(::Snowpack, state)

Retrieve the current snow density.
"""
snowdensity(::Snowpack, state) = state.ρsn

"""
    snowdepth(::Snowpack, ::SnowMassBalance, state)

Retrieve the current snow depth.
"""
snowdepth(snow::Snowpack, ::SnowMassBalance, state) = CryoGrid.thickness(snow, state)

"""
    snowfall(::SnowBC, state)

Retrieves the current snowfall flux for the given `SnowBC`. Defaults to `state.jw_snow`.
"""
snowfall(::SnowBC, state) = state.jw_snow

"""
    snowvariables(::Snowpack)

Default variables common to all snow schemes.
"""
snowvariables(::Snowpack) = (
    Diagnostic(:dsn, Scalar, u"m", domain=0..Inf),
    Diagnostic(:T_ub, Scalar, u"°C"),
)

### Default implmentations of CryoGrid methods for all Snowpack types ###

# for prescribed snow depth/density, the mass balance is given so we can skip computefluxes!
CryoGrid.computefluxes!(::Snowpack, ::SnowMassBalance, ssnow) = nothing

# Heat methods;
# thermal properties of snowpack
Heat.thermalproperties(snow::Snowpack) = snow.para.heat

# Hydrology methods;
Hydrology.hydraulicproperties(snow::Snowpack) = snow.para.water

# for snow, use constant hydraulic conductivity
Hydrology.hydraulicconductivity(snow::Snowpack, water::WaterBalance, θw, θwi, θsat) = Hydrology.kwsat(snow, water)

# max (fully saturated) water content
Hydrology.maxwater(::Snowpack, ::WaterBalance, state) = 1.0

# Default implementations of CryoGrid methods for Snowpack
CryoGrid.processes(snow::Snowpack) = Coupled(snow.mass, snow.water, snow.heat)

CryoGrid.thickness(::Snowpack, state, i::Integer=1) = abs(getscalar(state.Δz))

CryoGrid.midpoint(::Snowpack, state, i::Integer=1) = abs(getscalar(state.z) + getscalar(state.Δz)) / 2

CryoGrid.isactive(snow::Snowpack, state) = CryoGrid.thickness(snow, state) > threshold(snow)

# volumetric fractions for snowpack
function CryoGrid.volumetricfractions(::Snowpack, state, i)
    @inbounds let θwi = state.θwi[i],
        θw = state.θw[i],
        θa = 1.0 - θwi,
        θi = θwi - θw;
        return (θw, θi, θa)
    end
end

function CryoGrid.computediagnostic!(
    snow::Snowpack,
    procs::CoupledSnowWaterHeat,
    state,
)
    mass, water, heat = procs
    computediagnostic!(snow, mass, state)
    computediagnostic!(snow, water, state)
    computediagnostic!(snow, heat, state)
end

function CryoGrid.computediagnostic!(
    snow::Snowpack,
    mass::SnowMassBalance,
    state,
)
    # update snow density
    snowdensity!(snow, mass, state)
    # update snow depth;
    # by default, we just use the current layer thickness
    @setscalar state.dsn = getscalar(state.Δz)
end

# Special overrides for heat timestep control on snow layer

CryoGrid.timestep(::Snowpack, heat::HeatBalance{<:FreeWater,THeatOp,<:CryoGrid.CFL}, state) where {THeatOp} = error("CFL is not supported on snow layer")
function CryoGrid.timestep(snow::Snowpack, heat::HeatBalance{<:FreeWater,THeatOp,<:CryoGrid.MaxDelta}, state) where {THeatOp}
    Δx = Δ(state.grid)
    dtmax = Inf
    if getscalar(state.dsn) > snow.para.thresh
        @inbounds for i in eachindex(Δx)
            dtmax = min(dtmax, heat.dtlim(state.dH[i], state.H[i], state.t))
        end
        dtmax = isfinite(dtmax) && dtmax > 0 ? dtmax : Inf
    end
    return dtmax
end
