# Heat methods
"""
    thermalproperties(::SubSurface)

Returns the thermal properties for the given subsurface layer.
"""
function thermalproperties end
"""
    freezethaw!(sub::SubSurface, heat::HeatBalance, state)

Calculates freezing and thawing effects, including evaluation of the freeze curve.
In general, this function should compute at least the liquid/frozen water contents
and the corresponding heat capacity. Other variables such as temperature or enthalpy
may also need to be computed depending on the thermal scheme being implemented.
"""
function freezethaw! end

"""
    thermalconductivities(::SubSurface)

Get thermal conductivities for generic `SubSurface` layer.
"""
function thermalconductivities end

"""
    heatcapacities(::SubSurface)

Get heat capacities for generic `SubSurface` layer.
"""
function heatcapacities end

# Helper methods
"""
    enthalpy(T, C, L, θ) = T*C + L*θ

Discrete enthalpy function on temperature, heat capacity, specific latent heat of fusion, and liquid water content.
"""
@inline enthalpy(T, C, L, θ) = T*C + L*θ
"""
    enthalpyinv(H, C, L, θ) = (H - L*θ) / C

Discrete inverse enthalpy function given H, C, L, and θ.
"""
@inline enthalpyinv(H, C, L, θ) = (H - L*θ) / C
"""
    dHdT(T, C, L, ∂θw∂T, ch_w, ch_i) = C + ∂θw∂T*(L + T*(ch_w - ch_i))

Computes the apparent or "effective" heat capacity `∂H∂T` as a function of temperature, volumetric heat capacity,
latent heat of fusion, derivative of the freeze curve `∂θw∂T`, and the constituent heat capacities of water and ice.
"""
@inline dHdT(T, C, L, ∂θw∂T, ch_w, ch_i) = C + ∂θw∂T*(L + T*(ch_w - ch_i))
"""
    TemperatureProfile(pairs::Pair{<:Union{DistQuantity,Param},<:Union{TempQuantity,Param}}...)

Convenience constructor for `Numerics.Profile` which automatically converts temperature quantities.
"""
TemperatureProfile(pairs::Pair{<:Union{DistQuantity,Param},<:Union{TempQuantity,Param}}...) = Profile(map(p -> uconvert(u"m", p[1]) => uconvert(u"°C", p[2]), pairs))
