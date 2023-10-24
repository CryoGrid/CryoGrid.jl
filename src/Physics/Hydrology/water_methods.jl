"""
    waterdensity(sub::SubSurface)

Retrieves the density of water `ρw` from the given `SubSurface` layer. Default implementation assumes that
`WaterBalance` is provided as a field `water` on `sub`; this can of course, however, be overridden.
"""
waterdensity(sub::SubSurface) = sub.water.prop.ρw

"""
    hydraulicproperties(::SubSurface)

Retrieves the hydraulic properties from the given subsurface layer.
"""
hydraulicproperties(::SubSurface) = error("not implemented")
hydraulicproperties(sub::SubSurface, state) = hydraulicproperties(sub)
hydraulicproperties(sub::SubSurface, state, i) = hydraulicproperties(sub, state)

"""
    kwsat(::SubSurface, state, i)

Hydraulic conductivity at saturation.
"""
kwsat(sub::SubSurface) = hydraulicproperties(sub).kw_sat
kwsat(sub::SubSurface, state) = hydraulicproperties(sub, state).kw_sat
kwsat(sub::SubSurface, state, i) = hydraulicproperties(sub, state, i).kw_sat

"""
    interact_ET!(::Top, ::WaterBC, ::SubSurface, ::WaterBalance, state1, state2)
    interact_ET!(::SubSurface, ::WaterBalance, ::Bottom, ::WaterBC, state1, state2)
    interact_ET!(::SubSurface, ::WaterBalance, ::SubSurface, ::WaterBalance, state1, state2)

Specialized layer interaction for evapotranspirative processes. Default implementation does nothing.
Can be overridden by ET schemes and invoked in the relevant `interact!` implementation.
"""
interact_ET!(::Top, ::WaterBC, ::SubSurface, ::WaterBalance, state1, state2) = nothing
interact_ET!(::SubSurface, ::WaterBalance, ::Bottom, ::WaterBC, state1, state2) = nothing
interact_ET!(::SubSurface, ::WaterBalance, ::SubSurface, ::WaterBalance, state1, state2) = nothing

"""
    maxwater(sub::SubSurface, ::WaterBalance) 
    maxwater(sub::SubSurface, water::WaterBalance, state)
    maxwater(::SubSurface, ::WaterBalance, state, i)

Returns the maximum volumetric water content (saturation point) for grid cell `i`.
"""
maxwater(sub::SubSurface, ::WaterBalance) = error("maxwater not implemented for subsurface layer $sub")
maxwater(sub::SubSurface, water::WaterBalance, state) = maxwater(sub, water)
maxwater(sub::SubSurface, water::WaterBalance, state, i) = Utils.getscalar(maxwater(sub, water, state), i)

"""
    minwater(::SubSurface, water::WaterBalance)
    minwater(::SubSurface, water::WaterBalance, state, i)

Returns the minimum volumetric water content (typically field capacity for simplified schemes) for grid cell `i`. Defaults to zero.
"""
minwater(sub::SubSurface, ::WaterBalance) = sqrt(eps())
minwater(sub::SubSurface, water::WaterBalance, state) = minwater(sub, water)
minwater(sub::SubSurface, water::WaterBalance, state, i) = Utils.getscalar(minwater(sub, water, state), i)

"""
    watercontent(::SubSurface, state)
    watercontent(::SubSurface, state, i)

Returns the total water content `θwi` from the given subsurface layer and/or current state.
"""
watercontent(sub::SubSurface, state) = state.θwi
watercontent(sub::SubSurface, state, i) = Utils.getscalar(watercontent(sub, state), i)

"""
    hydraulicconductivity(water::WaterBalance, kw_sat, θw, θwi, θsat)

Computes the hydraulic conductivity for the given water balance configuration, current unfrozen
water content `θw`, total water/ice content `θwi`, and saturated (maximum) water content `θsat`.
"""
hydraulicconductivity(::WaterBalance, kw_sat, θw, θwi, θsat) = kw_sat*θw / θsat
