function heatcapacity(params::SoilParams, totalWater, liquidWater, mineral, organic)
    @unpack cw,co,cm,ca,ci = params.hc
    let air = 1.0 - totalWater - mineral - organic,
        ice = totalWater - liquidWater,
        water = liquidWater;
        water*cw + ice*ci + mineral*cm + organic*co + air*ca
    end
end
function thermalconductivity(params::SoilParams, totalWater, liquidWater, mineral, organic)
    @unpack kw,ko,km,ka,ki = params.tc
    let air = 1.0 - totalWater - mineral - organic,
        ice = totalWater - liquidWater,
        water = liquidWater;
        (water*kw^0.5 + ice*ki^0.5 + mineral*km^0.5 + organic*ko^0.5 + air*ka^0.5)^2
    end
end
""" Heat capacity for soil layer """
function heatcapacity!(soil::Soil, ::Heat, state)
    @. state.C = heatcapacity(soil.params, state.θw, state.θl, state.θm, state.θo)
end
""" Thermal conductivity for soil layer """
function thermalconductivity!(soil::Soil, ::Heat, state)
    @. state.kc = thermalconductivity(soil.params, state.θw, state.θl, state.θm, state.θo)
end

include("sfcc.jl")

""" Initial condition for heat conduction (all state configurations) on soil layer. """
function initialcondition!(soil::Soil, heat::Heat{U,<:SFCC}, state) where U
    interpolateprofile!(heat.profile, state)
    L = heat.params.L
    sfcc = freezecurve(heat)
    state.θl .= sfcc.f.(state.T, sfccparams(sfcc.f, soil, heat, state)...)
    @. state.C = heatcapacity(soil.params, state.θw, state.θl, state.θm, state.θo)
    @. state.H = enthalpy(state.T, state.C, L, state.θl)
end
""" Diagonstic step for heat conduction (all state configurations) on soil layer. """
function initialcondition!(soil::Soil, heat::Heat{(:Hₛ,:Hₗ),<:SFCC}, state)
    interpolateprofile!(heat.profile, state)
    L = heat.params.L
    sfcc = freezecurve(heat)
    state.θl .= sfcc.f.(state.T, sfccparams(sfcc.f, soil, heat, state)...)
    @. state.C = heatcapacity(soil.params, state.θw, state.θl, state.θm, state.θo)
    @. state.Hₛ = (state.T - 273.15)*state.C
    @. state.Hₗ = state.θl*L
end
