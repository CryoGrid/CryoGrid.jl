# === Thermal properties ===
# We use methods with optioinal index arguments `i` to allow for implementations both
# where these variables are treated as constants and as state variables.
# In the latter case, specializations should override only the index-free form
# and return a state vector instead of a scalar. The `getscalar` function will
# handle both the scalar and vector case!
"""
    mineral(soil::Soil, ::Heat, state)
    mineral(soil::Soil, ::Heat, state, i)

Retrieves the mineral content for the given layer at grid cell `i`, if provided.
Defaults to using the scalar mineral content defined on `soil`.
"""
@inline mineral(soil::Soil, ::Heat, state) = mineral(soil)
@inline mineral(soil::Soil, heat::Heat, state, i) = Utils.getscalar(mineral(soil, heat, state), i)
"""
    organic(soil::Soil, ::Heat, state)
    organic(soil::Soil, ::Heat, state, i)

Retrieves the organic content for the given layer at grid cell `i`, if provided.
Defaults to using the scalar organic content defined on `soil`.
"""
@inline organic(soil::Soil, ::Heat, state) = organic(soil)
@inline organic(soil::Soil, heat::Heat, state, i) = Utils.getscalar(organic(soil, heat, state), i)
"""
    porosity(soil::Soil, ::Heat, state)
    porosity(soil::Soil, ::Heat, state, i)

Retrieves the porosity for the given layer at grid cell `i`, if provided.
Defaults to using the scalar porosity defined on `soil`.
"""
@inline porosity(soil::Soil, ::Heat, state) = porosity(soil)
@inline porosity(soil::Soil, heat::Heat, state, i) = Utils.getscalar(porosity(soil, heat, state), i)
@inline heatcapacity(soil::Soil, heat::Heat, state, i) = heatcapacity(
    soil,
    totalwater(soil, heat, state, i),
    liquidwater(soil, heat, state, i),
    mineral(soil, heat, state, i),
    organic(soil, heat, state, i),
)
@inline function heatcapacity(soil::Soil, totalwater, liquidwater, mineral, organic)
    @unpack cw, co, cm, ca, ci = soil.hc
    let air = 1.0 - totalwater - mineral - organic,
        ice = totalwater - liquidwater,
        liq = liquidwater;
        liq*cw + ice*ci + mineral*cm + organic*co + air*ca
    end
end
@inline thermalconductivity(soil::Soil, heat::Heat, state, i) = thermalconductivity(
    soil,
    totalwater(soil, heat, state, i),
    liquidwater(soil, heat, state, i),
    mineral(soil, heat, state, i),
    organic(soil, heat, state, i),
)
@inline function thermalconductivity(soil::Soil, totalwater, liquidwater, mineral, organic)
    @unpack kw, ko, km, ka, ki = soil.tc
    let air = 1.0 - totalwater - mineral - organic,
        ice = totalwater - liquidwater,
        liq = liquidwater;
        (liq*kw^0.5 + ice*ki^0.5 + mineral*km^0.5 + organic*ko^0.5 + air*ka^0.5)^2
    end
end

# SFCC
include("sfcc.jl")

"""
Initial condition for heat conduction (all state configurations) on soil layer.
"""
function initialcondition!(soil::Soil, heat::Heat{<:SFCC}, state)
    initialcondition!(soil, state)
    L = heat.L
    sfcc = freezecurve(heat)
    state.θl .= sfcc.f.(state.T, sfccparams(sfcc.f, soil, heat, state)...)
    heatcapacity!(soil, heat, state)
    @. state.H = enthalpy(state.T, state.C, L, state.θl)
end
"""
Initial condition for heat conduction (all state configurations) with free water freeze curve on soil layer.
"""
function initialcondition!(soil::Soil, heat::Heat{FreeWater}, state)
    initialcondition!(soil, state)
    L = heat.L
    # initialize liquid water content based on temperature
    @inbounds for i in 1:length(state.T)
        θw = totalwater(soil, heat, state, i)
        state.θl[i] = ifelse(state.T[i] > 0.0, θw, 0.0)
        state.C[i] = heatcapacity(soil, heat, state, i)
        state.H[i] = enthalpy(state.T[i], state.C[i], L, state.θl[i])
    end
end
