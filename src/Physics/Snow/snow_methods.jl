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
ablation!(::Top, ::SnowBC, ::Snowpack, ::SnowMassBalance, stop, ssnow) = error("not implemented")

"""
    accumulate!(::Top, ::SnowBC, ::Snowpack, ::SnowMassBalance, stop, ssnow) 

Computes snow mass balance fluxes due to accumulation (e.g. snowfall).
"""
accumulate!(::Top, ::SnowBC, ::Snowpack, ::SnowMassBalance, stop, ssnow) = error("not implemented")

# Optional (w/ default implementations)

"""
    threshold(::Snowpack)

Retrieves the snow cover threshold for this `Snowpack` layer to become active.
"""
threshold(::Snowpack) = 0.01 # meters
"""
    swe(::Snowpack, ::SnowMassBalance, state)

Retrieve the current snow water equivalent of the snowpack.
"""
swe(::Snowpack, ::SnowMassBalance, state) = state.swe
# special implementation of swe for prescribed snow mass balance
swe(::Snowpack, smb::PrescribedSnowMassBalance, state) = smb.para.swe
swe(::Snowpack, smb::PrescribedSnowMassBalance{<:Forcing{u"m"}}, state) = smb.para.swe(state.t)

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

### Default implmentations of CryoGrid methods for Snowpack ###

# implement CryoGrid.Volume for prescribed vs. dynamic snow mass balance
CryoGrid.Volume(::Type{<:Snowpack{T,<:PrescribedSnowMassBalance}}) where {T} = CryoGrid.DiagnosticVolume()
CryoGrid.Volume(::Type{<:Snowpack{T,<:DynamicSnowMassBalance}}) where {T} = CryoGrid.PrognosticVolume()

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
