# Snow methods
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
swe(::Snowpack, smb::SnowMassBalance{<:PrescribedSnow}, state) = smb.para.swe
swe(::Snowpack, smb::SnowMassBalance{<:PrescribedSnow{<:Forcing{u"m"}}}, state) = smb.para.swe(state.t)

"""
    snowdensity(::Snowpack, ::SnowMassBalance, state)

Retrieve the current snow density.
"""
snowdensity(::Snowpack, ::SnowMassBalance, state) = state.ρsn

"""
    snowdensity!(::Snowpack, ::SnowMassBalance, state)

Computes the current snow density (if necessary) and stores the result in `state.ρsn`.
"""
snowdensity!(::Snowpack, ::SnowMassBalance, state) = error("not implemented")

"""
    snowdepth(::Snowpack, ::SnowMassBalance, state)

Retrieve the current snow depth.
"""
snowdepth(snow::Snowpack, ::SnowMassBalance, state) = CryoGrid.thickness(snow, state)

snowvariables(::Snowpack) = (
    Diagnostic(:dsn, Scalar, u"m", domain=0..Inf),
    Diagnostic(:T_ub, Scalar, u"°C"),
)
