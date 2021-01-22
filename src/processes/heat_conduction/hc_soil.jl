variables(soil::Soil, heat::Heat{UType"J"}) = (
    Var(:H, UFloat"J", OnGrid(Cells), Prognostic)
    Var(:T, UFloat"K", OnGrid(Cells)),
    Var(:dHdT, UFloat"J/(K*m^3)", OnGrid(Cells)),
    Var(:k, UFloat"W/(m*K)", OnGrid(Edges))
)

function enthalpyInv(heat::Heat{UType"J"}, H::UFloat"J", C::UFloat"J/(K*m^3)", totalWater)
    let ρ = heat.config.ρ,
        Lsl = heat.config.Lsl,
        L = ρ*Lsl, #[J/m^3]
        θ = max(1.0e-8, totalWater), #[Vol. fraction]
        # indicator variables for thawed and frozen states respectively
        I_t = @> H > L*θ float,
        I_f = @> H < 0.0 float;
        T = I_t*((H - L*θ) / C) + I_f*(H / C)
    end
end

"""
    enthalpy(T,water,hc)

Enthalpy at temperature T with the given water content and heat capacity.
"""
function enthalpy(heat::Heat{UType"J"}, T::UFloat"K", C::UFloat"J/(K*m^3)", liquidWater)
    let ρ = heat.config.ρ,
        Lsl = heat.config.Lsl,
        L = ρ*Lsl, #[J/m^3]
        θ = liquidWater; #[Vol. fraction]
        H = T*hc + θ*L
    end
end
