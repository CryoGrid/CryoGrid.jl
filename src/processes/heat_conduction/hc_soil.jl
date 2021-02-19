variables(soil::Soil, heat::Heat{UT"J"}) = (
    Var(:H, UFloat"J", OnGrid(Cells), Prognostic),
    Var(:T, UFloat"K", OnGrid(Cells)),
    Var(:C, UFloat"J/(K*m^3)", OnGrid(Cells)),
    Var(:k, UFloat"W/(m*K)", OnGrid(Edges))
)

function enthalpyInv(heat::Heat{UT"J"}, H::UFloat"J", C::UFloat"J/(K*m^3)", totalWater)
    let ρ = heat.config.ρ,
        Lsl = heat.config.Lsl,
        L = ρ*Lsl, #[J/m^3]
        θ = max(1.0e-8, totalWater), #[Vol. fraction]
        # indicator variables for thawed and frozen states respectively
        I_t = @> H > L*θ float;
        T = (H - I_t*L*θ) / C
    end
end

"""
    enthalpy(T,water,hc)

Enthalpy at temperature T with the given water content and heat capacity.
"""
function enthalpy(heat::Heat{UT"J"}, T::UFloat"K", C::UFloat"J/(K*m^3)", liquidWater)
    let ρ = heat.config.ρ,
        Lsl = heat.config.Lsl,
        L = ρ*Lsl, #[J/m^3]
        θ = liquidWater; #[Vol. fraction]
        H = T*hc + θ*L
    end
end

"""
Phase change with linear freeze curve. Assumes diagnostic liquid water variable. Should *not* be used with prognostic
water variable.
"""
function freezethaw(heat::Heat{UT"J",P}, totalWater, liquidWater) where {P<:HeatParams{LinearFC}}
    let ρ = heat.config.ρ,
        Lsl = heat.config.Lsl,
        L = ρ*Lsl, #[J/m^3]
        θ = max(1.0e-8, totalWater), #[Vol. fraction]
        Lθ = L*θ,
        I_t = @> H > Lθ float,
        I_c = @> H > 0.0 && H <= Lθ float;
        θ_l = I_c*(H/Lθ) + I_t
    end
end

# Note: @cryogrid macro not necessary here because we are in the same module; It's only here for consistency.
function diagnosticstep!(soil::Soil, heat::Heat{UT"J",P}, state) where {P<:HeatParams{LinearFC}}
    state.T .= enthalpyInv(heat, state.H, state.C, state.θ_w)
    state.θ_l .= freezethaw(heat, state.θ_w, state.θ_l)
end

function prognosticstep!(soil::Soil, heat::Heat{UT"J"}, state)
    δₓ = Δ(state.grids.T)
    # upper boundary
    state.dH[1] += let T₂=state.T[2], T₁=state.T[1], k=state.k[2], δ=δₓ[1];
        k*(T₂-T₁)/δ
    end
    # diffusion on non-boundary cells
    ∇²(state.T, δₓ, state.k, state.dH[2:end-1])
    # lower boundary
    state.dH[end] += let T₂=state.T[end], T₁=state.T[end-1], k=state.k[end-1], δ=δₓ[end];
        -k*(T₂-T₁)/δ
    end
    return nothing # ensure no allocation
end

"""
Top interaction, constant temperature Dirichlet boundary condition.
"""
function interact!(top::Top, c::ConstantAirTemp, soil::Soil, heat::Heat{UT"J"}, stop, ssoil)
    δₓ = ssoil.grids.T[1] - ssoil.grids.k[1] # distance to surface
    ssoil.dH[1] += let Tair=c.value, Tsoil=sstate.T[1], k=sstate.k[1];
        ∂H - k*(Tsoil-Tair)/δₓ
    end
    return nothing # ensure no allocation
end

function interact!(soil::Soil, heat::Heat{UT"J"}, bottom::Bottom, Qgeo::GeothermalHeatFlux, ssoil, sbot)
    ssoil.dH[end] += Qgeo.value
    return nothing # ensure no allocation
end
