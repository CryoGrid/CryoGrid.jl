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
    compaction!(::Top, ::SnowBC, ::Snowpack, ::SnowMassBalance, stop, ssnow)

Computes snow density changes due to compaction, if defined.
"""
compaction!(::Top, ::SnowBC, ::Snowpack, ::SnowMassBalance, stop, ssnow) = error("not implemented")

"""
    accumulation!(::Top, ::SnowBC, ::Snowpack, ::SnowMassBalance, stop, ssnow) 

Computes snow mass balance fluxes due to accumulation (e.g. snowfall).
"""
accumulation!(::Top, ::SnowBC, ::Snowpack, ::SnowMassBalance, stop, ssnow) = error("not implemented")

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
