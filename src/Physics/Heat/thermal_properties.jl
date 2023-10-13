# Define thermal properties type and constructor

"""
    ThermalProperties

Material thermal properties, e.g. conductivity and heat capacity. By default,
this includes the thermal properties of water, ice, and air. This can be extended
by passing additional properties into the constructor.
"""
Utils.@properties ThermalProperties(
    kh_w = 0.57u"J/s/m/K", # thermal conductivity of water [Hillel (1982)]
    kh_i = 2.2u"J/s/m/K", # thermal conductivity of ice [Hillel (1982)]
    kh_a = 0.025u"J/s/m/K", # thermal conductivity of air [Hillel (1982)]
    ch_w = 4.2e6u"J/K/m^3", # heat capacity of water
    ch_i = 1.9e6u"J/K/m^3", # heat capacity of ice
    ch_a = 0.00125e6u"J/K/m^3", # heat capacity of air
)

# Thermal conductivity

"""
    thermalconductivities(sub::SubSurface)

Get thermal conductivities for generic `SubSurface` layer.
"""
function thermalconductivities(sub::SubSurface)
    @unpack kh_w, kh_i, kh_a = thermalproperties(sub)
    return (kh_w, kh_i, kh_a)
end
thermalconductivities(sub::SubSurface, state) = thermalconductivities(sub)
thermalconductivities(sub::SubSurface, state, i) = thermalconductivities(sub, state)

"""
    quadratic_parallel_conductivity(ks, θs)

The "quadratic parallel" thermal conductivity formula as defined by Cosenza et al. 2003:

```math
k = [\\sum_{i=1}^N θᵢ\\sqrt{kᵢ}]^2
```

Cosenza, P., Guérin, R., and Tabbagh, A.: Relationship between thermal
conductivity and water content of soils using numerical modelling,
European Journal of Soil Science, 54, 581–588,
https://doi.org/10.1046/j.1365-2389.2003.00539.x, 2003.
"""
quadratic_parallel_conductivity(ks, θs) = sum(map(*, map(sqrt, ks), θs))^2

"""
    geometric_conductivity(ks::NTuple{N}, θs::NTuple{N}) where {N}

Geometric mean of constituent thermal conductivities according to Woodside and Messmer (1961).

Woodside, W. & Messmer, J.H. 1961. Thermal conductivity of porous media. I. Unconsolidated sands. Journal of Applied Physics, 32, 1688–1699. 
"""
geometric_conductivity(ks, θs) = prod(map((k,θ) -> k^θ, ks, θs))

"""
    thermalconductivity(sub::SubSurface, heat::HeatBalance, state, i)

Computes the thermal conductivity as a squared weighted sum over constituent conductivities with volumetric fractions `θfracs`.
"""
function thermalconductivity(sub::SubSurface, heat::HeatBalance, state, i)
    θs = volumetricfractions(sub, state, i)
    ks = thermalconductivities(sub, state, i)
    f = thermal_conductivity_function(heat.op)
    return f(ks, θs)
end

"""
    thermalconductivity!(sub::SubSurface, heat::HeatBalance, state)

Computes the thermal conductivity for the given layer from the current state and stores the result in-place in the state variable `k`.
"""
function thermalconductivity!(sub::SubSurface, heat::HeatBalance, state)
    @inbounds for i in 1:length(state.T)
        state.kc[i] = thermalconductivity(sub, heat, state, i)
        if i > 1
            Δk₁ = CryoGrid.thickness(sub, state, i-1)
            Δk₂ = CryoGrid.thickness(sub, state, i)
            state.k[i] = Numerics.harmonicmean(state.kc[i-1], state.kc[i], Δk₁, Δk₂)
        end
    end
    # thermal conductivity at boundaries
    # assumes boundary conductivities = cell conductivities
    @inbounds state.k[1] = state.kc[1]
    @inbounds state.k[end] = state.kc[end]
end

# Heat capacity

"""
    heatcapacities(::SubSurface)

Get heat capacities for generic `SubSurface` layer.
"""
function heatcapacities(sub::SubSurface)
    @unpack ch_w, ch_i, ch_a = thermalproperties(sub)
    return ch_w, ch_i, ch_a
end
heatcapacities(sub::SubSurface, state) = heatcapacities(sub)
heatcapacities(sub::SubSurface, state, i) = heatcapacities(sub, state)

"""
    weighted_average_heatcapacity(cs, θs)

Represents a simple composite heat capacity that is the sum of each constituent heat capacity weighted
by the volumetric fraction.
"""
weighted_average_heatcapacity(cs, θs) = sum(map(*, cs, θs))

"""
    heatcapacity(sub::SubSurface, heat::HeatBalance, state, i)

Computes the heat capacity as a weighted average over constituent capacities with volumetric fractions `θfracs`.
"""
function heatcapacity(sub::SubSurface, heat::HeatBalance, state, i)
    θs = volumetricfractions(sub, state, i)
    cs = heatcapacities(sub, state, i)
    f = heat_capacity_function(heat.op)
    return f(cs, θs)
end

"""
    heatcapacity!(sub::SubSurface, heat::HeatBalance, state)

Computes the heat capacity for the given layer from the current state and stores the result in-place in the state variable `C`.
"""
function heatcapacity!(sub::SubSurface, heat::HeatBalance, state)
    @inbounds for i in 1:length(state.T)
        state.C[i] = heatcapacity(sub, heat, state, i)
    end
end