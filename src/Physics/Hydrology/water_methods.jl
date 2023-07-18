"""
    hydraulicproperties(::SubSurface)

Retrieves the hydraulic properties from the given subsurface layer. Default implementation
simply returns the default configuration of `HydraulicProperties`.
"""
hydraulicproperties(::SubSurface) = error("not implemented")

"""
    waterdensity(sub::SubSurface)

Retrieves the density of water `ρw` from the given `SubSurface` layer. Default implementation assumes that
`WaterBalance` is provided as a field `water` on `sub`; this can of course, however, be overridden.
"""
waterdensity(sub::SubSurface) = sub.water.prop.ρw

"""
    kwsat(::SubSurface, ::WaterBalance)

Hydraulic conductivity at saturation.
"""
kwsat(sub::SubSurface, ::WaterBalance) = hydraulicproperties(sub).kw_sat

"""
    interact_ET!(::SubSurface, ::WaterBalance, ::SubSurface, ::WaterBalance, state1, state2)

Specialized subsurface layer interaction for evapotranspirative processes. Default implementation does nothing.
"""
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
minwater(sub::SubSurface, ::WaterBalance) = hydraulicproperties(sub).fieldcapacity
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
    watercontent!(::SubSurface, ::WaterBalance, state)

Computes the volumetric water content from current saturation or pressure state.
"""
@inline function watercontent!(sub::SubSurface, water::WaterBalance, state)
    @inbounds for i in eachindex(state.sat)
        state.θsat[i] = maxwater(sub, water, state, i)
        state.θwi[i] = state.sat[i]*state.θsat[i]
    end
end

"""
    hydraulicconductivity(sub::SubSurface, water::WaterBalance, θw, θwi, θsat)

Computes the hydraulic conductivity for the given layer and water balance configuration, current unfrozen
water content `θw`, total water/ice content `θwi`, and saturated (maximum) water content `θsat`.
"""
hydraulicconductivity(sub::SubSurface, water::WaterBalance, θw, θwi, θsat) = kwsat(sub, water)*θw / θsat

"""
    hydraulicconductivity!(sub::SubSurface, water::WaterBalance, state)

Computes hydraulic conductivities for the given subsurface layer and water balance scheme.
"""
@inline function hydraulicconductivity!(sub::SubSurface, water::WaterBalance, state)
    @inbounds for i in eachindex(state.kwc)
        let θsat = Hydrology.maxwater(sub, water, state, i),
            θw = state.θw[i],
            θwi = state.θwi[i];
            state.kwc[i] = hydraulicconductivity(sub, water, θw, θwi, θsat)
            if i > 1
                state.kw[i] = min(state.kwc[i-1], state.kwc[i])
            end
        end
    end
    # set hydraulic conductivity at boundaries equal to cell conductivities
    state.kw[1] = state.kwc[1]
    state.kw[end] = state.kwc[end]
end
