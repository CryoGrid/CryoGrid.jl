# === Thermal properties ===
# We use methods with optional index arguments `i` to allow for implementations both
# where these variables are treated as constants and as state variables.
# In the latter case, specializations should override only the index-free form
# and return a state vector instead of a scalar. The `getscalar` function will
# handle both the scalar and vector case!
"""
    mineral(soil::Soil, state, i)

Retrieves the mineral content for the given layer at grid cell `i`, if provided.
Defaults to using the scalar mineral content defined on `soil`.
"""
@inline mineral(soil::Soil, state, i) = Utils.getscalar(mineral(soil, state), i)
"""
    organic(soil::Soil, state, i)

Retrieves the organic content for the given layer at grid cell `i`, if provided.
Defaults to using the scalar organic content defined on `soil`.
"""
@inline organic(soil::Soil, state, i) = Utils.getscalar(organic(soil, state), i)
"""
    porosity(soil::Soil, state, i)

Retrieves the porosity for the given layer at grid cell `i`, if provided.
Defaults to using the scalar porosity defined on `soil`.
"""
@inline porosity(soil::Soil, state, i) = Utils.getscalar(porosity(soil, state), i)
# Functions for retrieving constituents and volumetric fractions
@inline thermalconductivities(soil::Soil, heat::Heat) = (heat.prop.kw, heat.prop.ki, heat.prop.ka, soil.prop.km, soil.prop.ko)
@inline heatcapacities(soil::Soil, heat::Heat) = (heat.prop.cw, heat.prop.ci, heat.prop.ca, soil.prop.cm, soil.prop.co)
@inline function volumetricfractions(soil::Soil, heat::Heat, state, i)
    return let θw = totalwater(soil, state, i),
        θl = liquidwater(soil, heat, state, i),
        θm = mineral(soil, state, i),
        θo = organic(soil, state, i),
        θa = 1.0 - θw - θm - θo,
        θi = θw - θl;
        (θl, θi, θa, θm, θo)
    end
end
@inline function heatcapacity(soil::Soil, heat::Heat, θw, θl, θm, θo)
    let θa = 1.0 - θw - θm - θo,
        θi = θw - θl;
        return heatcapacity(heatcapacities(soil, heat), (θl, θi, θa, θm, θo))
    end
end
@inline function thermalconductivity(soil::Soil, heat::Heat, θw, θl, θm, θo)
    let θa = 1.0 - θw - θm - θo,
        θi = θw - θl;
        return thermalconductivity(thermalconductivities(soil, heat), (θl, θi, θa, θm, θo))
    end
end

# SFCC
include("sfcc.jl")

"""
Initial condition for heat conduction (all state configurations) on soil layer w/ SFCC.
"""
function initialcondition!(soil::Soil, heat::Heat{<:SFCC}, state)
    initialcondition!(soil, heat, freezecurve(heat), state)
end
"""
Initial condition for heat conduction (all state configurations) on soil layer w/ free water freeze curve.
"""
function initialcondition!(soil::Soil, heat::Heat{FreeWater}, state)
    L = heat.L
    # initialize liquid water content based on temperature
    @inbounds for i in 1:length(state.T)
        θw = totalwater(soil, state, i)
        state.θl[i] = ifelse(state.T[i] > 0.0, θw, 0.0)
        state.C[i] = heatcapacity(soil, heat, state, i)
        state.H[i] = enthalpy(state.T[i], state.C[i], L, state.θl[i])
    end
end
