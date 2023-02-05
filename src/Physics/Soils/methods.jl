"""
    mineral(soil::Soil, state, i)

Retrieves the mineral content for the given layer at grid cell `i`, if provided.
Defaults to using the scalar mineral content defined on `soil`.
"""
mineral(soil::Soil, state, i) = Utils.getscalar(mineral(soil, state), i)
mineral(::Soil, state) = state.θm

"""
    organic(soil::Soil, state, i)

Retrieves the organic content for the given layer at grid cell `i`, if provided.
Defaults to using the scalar organic content defined on `soil`.
"""
organic(soil::Soil, state, i) = Utils.getscalar(organic(soil, state), i)
organic(::Soil, state) = state.θo

"""
    porosity(soil::Soil, state, i)

Retrieves the porosity for the given layer at grid cell `i`, if provided.
Defaults to using the scalar porosity defined on `soil`.
"""
porosity(soil::Soil, state, i) = Utils.getscalar(porosity(soil, state), i)
porosity(::Soil, state) = state.θp
