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
snowdensity(::Snowpack, smb::SnowMassBalance{<:PrescribedSnow}, state) = smb.para.ρsn
snowdensity(::Snowpack, smb::SnowMassBalance{<:PrescribedSnow{Tswe,<:Forcing{u"kg/m^3"}}}, state) where {Tswe} = smb.para.ρsn(state.t)

"""
    accumulation(snow::SnowMassBalance{<:DynamicSnow})

Get the snow accumulation scheme from the given `SnowMassBalance` parameterization.
"""
accumulation(snow::SnowMassBalance{<:DynamicSnow}) = snow.para.accumulation

"""
    ablation(snow::SnowMassBalance{<:DynamicSnow})

Get the snow ablation scheme from the given `SnowMassBalance` parameterization.
"""
ablation(snow::SnowMassBalance{<:DynamicSnow}) = snow.para.ablation

"""
    density(snow::SnowMassBalance{<:DynamicSnow})

Get the snow density scheme from the given `SnowMassBalance` parameterization.
"""
density(snow::SnowMassBalance{<:DynamicSnow}) = snow.para.density

# Default implementations of CryoGrid methods for Snowpack
CryoGrid.basevariables(::Snowpack, ::SnowMassBalance) = (
    Diagnostic(:dsn, Scalar, u"m", domain=0..Inf),
    Diagnostic(:T_ub, Scalar, u"°C"),
)
# for prescribed snow depth/density, the mass balance is given so we do not need to do anything here
CryoGrid.prognosticstep!(::Snowpack, ::SnowMassBalance{<:PrescribedSnow}, ssnow) = nothing
# thermal properties snowpack
Heat.thermalproperties(snow::Snowpack) = snow.prop.heat
# volumetric fractions for snowpack
@inline function CryoGrid.volumetricfractions(::Snowpack, state, i)
    @inbounds let θwi = state.θwi[i],
        θw = state.θw[i],
        θa = 1.0 - θwi,
        θi = θwi - θw;
        return (θw, θi, θa)
    end
end
