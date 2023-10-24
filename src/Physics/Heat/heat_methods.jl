# Heat methods
"""
    thermalproperties(::SubSurface)

Returns the thermal properties for the given subsurface layer.
"""
thermalproperties(::SubSurface) = error("not implemented")
thermalproperties(sub::SubSurface, state) = thermalproperties(sub)
thermalproperties(sub::SubSurface, state, i) = thermalproperties(sub, state)

"""
    thermalconductivity(::SubSurface, state, i)

Computes the thermal conductivity for the given `SubSurface` layer at grid cell `i`.
"""
thermalconductivity(::SubSurface, state, i) = error("not implemented")

"""
    heatcapacity(::SubSurface, state, i)

Computes the heat capacity for the given `SubSurface` layer at grid cell `i`.
"""
heatcapacity(::SubSurface, state, i) = error("not implemented")

"""
    freezethaw!(sub::SubSurface, heat::HeatBalance, state)

Calculates freezing and thawing effects, including evaluation of the freeze curve.
In general, this function should compute at least the liquid/frozen water contents
and the corresponding heat capacity. Other variables such as temperature or enthalpy
may also need to be computed depending on the thermal scheme being implemented.
"""
freezethaw!(::SubSurface, ::HeatBalance, state) = error("not implemented")

# Helper methods
"""
    enthalpy(T, C, L, θ) = T*C + L*θ

Discrete enthalpy function on temperature, heat capacity, specific latent heat of fusion, and liquid water content.
"""
enthalpy(T, C, L, θ) = T*C + L*θ

"""
    enthalpyinv(H, C, L, θ) = (H - L*θ) / C

Discrete inverse enthalpy function given H, C, L, and θ.
"""
enthalpyinv(H, C, L, θ) = (H - L*θ) / C

"""
    dHdT(T, C, L, ∂θw∂T, ch_w, ch_i) = C + ∂θw∂T*(L + T*(ch_w - ch_i))

Computes the apparent or "effective" heat capacity `∂H∂T` as a function of temperature, volumetric heat capacity,
latent heat of fusion, derivative of the freeze curve `∂θw∂T`, and the constituent heat capacities of water and ice.
"""
dHdT(T, C, L, ∂θw∂T, ch_w, ch_i) = C + ∂θw∂T*(L + T*(ch_w - ch_i))

"""
    TemperatureProfile(pairs::Pair{<:Union{DistQuantity,Param},<:Union{TempQuantity,Param}}...)

Convenience constructor for `Numerics.Profile` which automatically converts temperature quantities.
"""
TemperatureProfile(pairs::Pair{<:Union{DistQuantity,Param},<:Union{TempQuantity,Param}}...) = Profile(map(p -> uconvert(u"m", p[1]) => uconvert(u"°C", p[2]), pairs))

"""
    thermal_conductivity_function(op::HeatOperator)

Retreives the thermal conductivity function for this heat operator.
"""
thermal_conductivity_function(op::HeatOperator) = op.cond

"""
    heat_capacity_function(op::HeatOperator)

Retreives the volumetric heat capacity function for this heat operator.
"""
heat_capacity_function(op::HeatOperator) = op.hc
