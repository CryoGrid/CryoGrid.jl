# Composition
"""
    volumetricfractions(::Layer, ::SubSurfaceProcess, state)
    volumetricfractions(::Layer, ::SubSurfaceProcess, state, i)

Get the volumetric fractions of each constituent in the volume (at grid cell `i`, if specificed).
"""
volumetricfractions(::Layer, ::SubSurfaceProcess, state) = ()
volumetricfractions(::Layer, ::SubSurfaceProcess, state, i) = ()
"""
    waterice(::SubSurface, state)
    waterice(sub::SubSurface, ::SubSurfaceProcess, state)
    waterice(sub::SubSurface, ::SubSurfaceProcess, state, i)

Retrieves the total water content (water + ice) for the given layer at grid cell `i`, if specified.
Defaults to retrieving the state variable `θwi` (assuming it exists).
"""
@inline waterice(::SubSurface, state) = state.θwi
@inline waterice(sub::SubSurface, proc::SubSurfaceProcess, state) = waterice(sub, proc, state)
@inline waterice(sub::SubSurface, proc::SubSurfaceProcess, state, i) = Utils.getscalar(waterice(sub, proc, state), i)
