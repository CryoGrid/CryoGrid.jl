"""
    heatcapacities(::SubSurface, heat::HeatBalance)

Get heat capacities for generic `SubSurface` layer.
"""
@inline function heatcapacities(sub::SubSurface)
    @unpack hc_w, hc_i, hc_a = thermalproperties(sub)
    return hc_w, hc_i, hc_a
end
"""
    weighted_average_heatcapacity(cs, θs)

Represents a simple composite heat capacity that is the sum of each constituent heat capacity weighted
by the volumetric fraction.
"""
weighted_average_heatcapacity(cs, θs) = sum(map(*, cs, θs))
"""
    heatcapacity(sub::SubSurface, heat::HeatBalance, θfracs...)

Computes the heat capacity as a weighted average over constituent capacities with volumetric fractions `θfracs`.
"""
@inline function heatcapacity(sub::SubSurface, heat::HeatBalance, θfracs...)
    cs = heatcapacities(sub)
    f = heatcapacity(operator(heat))
    return f(cs, θfracs)
end
"""
    heatcapacity!(sub::SubSurface, heat::HeatBalance, state)

Computes the heat capacity for the given layer from the current state and stores the result in-place in the state variable `C`.
"""
@inline function heatcapacity!(sub::SubSurface, heat::HeatBalance, state)
    @inbounds for i in 1:length(state.T)
        θfracs = volumetricfractions(sub, state, i)
        state.C[i] = heatcapacity(sub, heat, θfracs...)
    end
end
