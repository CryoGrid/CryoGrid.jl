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
