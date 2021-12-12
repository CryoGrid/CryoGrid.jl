@inline function heatcapacity(soil::Soil, totalwater, liquidwater, mineral, organic)
    @unpack cw,co,cm,ca,ci = soil.hc
    let air = 1.0 - totalwater - mineral - organic,
        ice = totalwater - liquidwater,
        liq = liquidwater;
        liq*cw + ice*ci + mineral*cm + organic*co + air*ca
    end
end
@inline function thermalconductivity(soil::Soil, totalwater, liquidwater, mineral, organic)
    @unpack kw,ko,km,ka,ki = soil.tc
    let air = 1.0 - totalwater - mineral - organic,
        ice = totalwater - liquidwater,
        liq = liquidwater;
        (liq*kw^0.5 + ice*ki^0.5 + mineral*km^0.5 + organic*ko^0.5 + air*ka^0.5)^2
    end
end
""" Heat capacity for soil layer """
@inline function heatcapacity!(soil::Soil, ::Heat, state)
    @inbounds @. state.C = heatcapacity(soil, state.θw, state.θl, state.θm, state.θo)
end
""" Thermal conductivity for soil layer """
@inline function thermalconductivity!(soil::Soil, ::Heat, state)
    @inbounds @. state.kc = thermalconductivity(soil, state.θw, state.θl, state.θm, state.θo)
end

include("sfcc.jl")

""" Initial condition for heat conduction (all state configurations) on soil layer. """
function initialcondition!(soil::Soil, heat::Heat{<:SFCC}, state)
    initialcondition!(soil, state)
    L = heat.L
    sfcc = freezecurve(heat)
    state.θl .= sfcc.f.(state.T, sfccparams(sfcc.f, soil, heat, state)...)
    heatcapacity!(soil, heat, state)
    @. state.H = enthalpy(state.T, state.C, L, state.θl)
end
