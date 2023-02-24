"""
    kwsat(::SubSurface, ::WaterBalance)

Hydraulic conductivity at saturation.
"""
function kwsat end
"""
    maxwater(::SubSurface, ::WaterBalance, state, i)

Returns the maximum volumetric water content (saturation point) for grid cell `i`. Defaults to `1`.
"""
function maxwater end
"""
    minwater(::SubSurface, water::WaterBalance)

Returns the minimum volumetric water content (typically field capacity for simplified schemes) for grid cell `i`. Defaults to zero.
"""
function minwater end
"""
    watercontent(::SubSurface, state)
    watercontent(::SubSurface, ::Process, state)
    watercontent(::SubSurface, state, i)

Returns the total water content `θwi` from the given subsurface layer and/or current state.
"""
function watercontent end
"""
    watercontent!(::SubSurface, ::WaterBalance, state)

Computes the volumetric water content from current saturation or pressure state.
"""
function watercontent! end
"""
    hydraulicconductivity!(sub::SubSurface, water::WaterBalance, state)

Computes hydraulic conductivities for the given subsurface layer and water balance scheme.
"""
function hydraulicconductivity! end
"""
    wateradvection!(sub::SubSurface, water::WaterBalance, state)

Computes the advective component of water fluxes due to gravity and stores the result in `state.jw`.
"""
function wateradvection! end
"""
    waterdiffusion!(::SubSurface, ::WaterBalance, state)

Computes diffusive fluxes for water balance, if defined.
"""
function waterdiffusion! end
"""
    waterprognostic!(::SubSurface, ::WaterBalance, state)

Computes the prognostic time derivative for the water balance, usually based on `∂θwi∂t`.
Implementation depends on which water flow scheme is being used.
"""
function waterprognostic! end
