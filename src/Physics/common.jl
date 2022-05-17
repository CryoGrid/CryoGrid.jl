# Composition
"""
    volumetricfractions(::Layer, ::Process, state)
    volumetricfractions(::Layer, ::Process, state, i)

Get the volumetric fractions of each constituent in the volume (at grid cell `i`, if specificed).
"""
volumetricfractions(::Layer, ::Process, state) = ()
volumetricfractions(::Layer, ::Process, state, i) = ()
"""
    waterice(sub::SubSurface, state)
    waterice(sub::SubSurface, state, i)

Retrieves the total water content (water + ice) for the given layer at grid cell `i`, if specified.
Defaults to retrieving the state variable `θwi` (assuming it exists).
"""
@inline waterice(::SubSurface, state) = state.θwi
@inline waterice(sub::SubSurface, state, i) = Utils.getscalar(waterice(sub, state), i)
