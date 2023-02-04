"""
    thermalconductivities(::SubSurface, heat::HeatBalance)

Get thermal conductivities for generic `SubSurface` layer.
"""
@inline function thermalconductivities(sub::SubSurface)
    @unpack kh_w, kh_i, kh_a = thermalproperties(sub)
    return (kh_w, kh_i, kh_a)
end
"""
    thermalconductivity(op::HeatOperator)

Returns the thermal conductivity function for the given `HeatOperator`.
"""
thermalconductivity(op::HeatOperator) = error("not implemented for $(typeof(op))")
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
    thermalconductivity(sub::SubSurface, heat::HeatBalance, θfracs...)

Computes the thermal conductivity as a squared weighted sum over constituent conductivities with volumetric fractions `θfracs`.
"""
@inline function thermalconductivity(sub::SubSurface, heat::HeatBalance, θfracs...)
    ks = thermalconductivities(sub)
    f = thermalconductivity(operator(heat))
    return f(ks, θfracs)
end
"""
    thermalconductivity!(sub::SubSurface, heat::HeatBalance, state)

Computes the thermal conductivity for the given layer from the current state and stores the result in-place in the state variable `C`.
"""
@inline function thermalconductivity!(sub::SubSurface, heat::HeatBalance, state)
    @inbounds for i in 1:length(state.T)
        θfracs = volumetricfractions(sub, state, i)
        state.kc[i] = thermalconductivity(sub, heat, θfracs...)
    end
    # thermal conductivity at boundaries
    # assumes boundary conductivities = cell conductivities
    @inbounds state.k[1] = state.kc[1]
    @inbounds state.k[end] = state.kc[end]
    # Harmonic mean of inner conductivities
    @inbounds let k = (@view state.k[2:end-1]),
        Δk = Δ(state.grids.k);
        Numerics.harmonicmean!(k, state.kc, Δk)
    end
end