# Heat methods
"""
    thermalproperties(::SubSurface)

Returns the thermal properties for the given subsurface layer.
"""
thermalproperties(::SubSurface) = error("not implemented")
thermalproperties(sub::SubSurface, state) = thermalproperties(sub)
thermalproperties(sub::SubSurface, state, i) = thermalproperties(sub, state)

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
    thermalconductivity(sub::SubSurface, heat::HeatBalance, state, i)

Computes the thermal conductivity as a squared weighted sum over constituent conductivities with volumetric fractions `θfracs`.
"""
function thermalconductivity(sub::SubSurface, state, i)
    θs = volumetricfractions(sub, state, i)
    ks = thermalconductivities(sub, state, i)
    f = thermalproperties(sub).thermcond
    return f(ks, θs)
end

"""
    heatcapacity(sub::SubSurface, state, i)

Computes the heat capacity as a weighted average over constituent capacities with volumetric fractions `θfracs`.
"""
function heatcapacity(sub::SubSurface, state, i)
    θs = volumetricfractions(sub, state, i)
    cs = heatcapacities(sub, state, i)
    f = thermalproperties(sub).heatcap
    return f(cs, θs)
end

"""
    freezecurve(sub::SubSurface)

Returns the soil freezing characteristic `FreezeCurve` for the given subsurface layer.
Defautls to `FreeWater`.
"""
freezecurve(sub::SubSurface) = FreeWater()

"""
    freezethaw!(::FreezeCurve, ::SubSurface, ::Process, state)

Calculates freezing and thawing effects, including evaluation of the freeze curve.
In general, this function should compute at least the liquid/frozen water contents
and the corresponding heat capacity. Other variables such as temperature or enthalpy
may also need to be computed depending on the thermal scheme being implemented.
"""
freezethaw!(::FreezeCurve, ::SubSurface, ::Process, state) = error("not implemented")

"""
    enthalpyinv([::FreezeCurve], sub::SubSurface, heat::HeatBalance, state)

Evaluates the inverse enthalpy function (H -> T) on the current state.
"""
enthalpyinv(sub::SubSurface, heat::HeatBalance, state) = enthalypinv(freezecurve(sub), sub, heat, state)

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
TemperatureProfile(pairs::Pair{<:Union{DistQuantity,Param}}...) = Profile(map(p -> uconvert(u"m", p[1]) => p[2], pairs))
