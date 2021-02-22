variables(soil::Soil, heat::Heat{UT"J"}) = (
    Prognostic(:H, UFloat"J", OnGrid(Cells)),
    Diagnostic(:T, UFloat"K", OnGrid(Cells)),
    Diagnostic(:C, UFloat"J/(K*m^3)", OnGrid(Cells)),
    Diagnostic(:k, UFloat"W/(m*K)", OnGrid(Edges)),
    Diagnostic(:θw_k, UFloat"W/(m*K)", OnGrid(Edges)),
    Diagnostic(:θl_k, UFloat"W/(m*K)", OnGrid(Edges)),
    Diagnostic(:θm_k, UFloat"W/(m*K)", OnGrid(Edges)),
    Diagnostic(:θo_k, UFloat"W/(m*K)", OnGrid(Edges)),
)

function enthalpyInv(H::UFloat"J", C::UFloat"J/(K*m^3)", totalWater, L::UFloat"J/m^3")
    let θ = max(1.0e-8, totalWater), #[Vol. fraction]
        # indicator variables for thawed and frozen states respectively
        I_t = (H > L*θ) |> float;
        T = (H - I_t*L*θ) / C
    end
end

"""
    enthalpy(T,water,hc)

Enthalpy at temperature T with the given water content and heat capacity.
"""
function enthalpy(T::UFloat"K", C::UFloat"J/(K*m^3)", liquidWater, L::UFloat"J/m^3")
    let θ = liquidWater; #[Vol. fraction]
        H = T*C + θ*L
    end
end

"""
Phase change with linear freeze curve. Assumes diagnostic liquid water variable. Should *not* be used with prognostic
water variable.
"""
function freezethaw(H::UFloat"J", totalWater, L::UFloat"J/m^3")
    let θ = max(1.0e-8, totalWater), #[Vol. fraction]
        Lθ = L*θ,
        I_t = (H > Lθ) |> float,
        I_c = (H > 0.0 && H <= Lθ) |> float;
        liquidfraction = I_c*(H/Lθ) + I_t
    end
end

export enthalpy, enthalpyInv, freezethaw

function initialcondition!(soil::Soil, heat::Heat{UT"J"}, state)
    interpolateprofile!(heat.profile, state)
    L = ρ(heat)*Lsl(heat)
    @. state.C = heatCapacity(soil.hcparams, state.θw, state.θl, state.θm, state.θo)
    @. state.H = enthalpy(state.T, state.C, state.θl, L)
    # k lies on the boundary grid, so we have to regrid the soil properties
    regrid!(state.θw_k, state.θw, state.grids.θw, state.grids.k)
    regrid!(state.θl_k, state.θl, state.grids.θl, state.grids.k)
    regrid!(state.θm_k, state.θm, state.grids.θm, state.grids.k)
    regrid!(state.θo_k, state.θo, state.grids.θo, state.grids.k)
    @. state.k = thermalConductivity(soil.tcparams, state.θw_k, state.θl_k, state.θm_k, state.θo_k)
end

function diagnosticstep!(soil::Soil, heat::Heat{UT"J",P}, state) where {P<:HeatParams{LinearFC}}
    let ρ = heat.params.ρ,
        Lsl = heat.params.Lsl,
        L = ρ*Lsl; #[J/m^3];
        @. state.T = enthalpyInv(state.H, state.C, state.θw, L)
        @. state.θl = freezethaw(state.H, state.θw, L)
        @. state.θl = state.θw*state.θl
    end
    @. state.C = heatCapacity(soil.hcparams, state.θw, state.θl, state.θm, state.θo)
    regrid!(state.θl_k, state.θl, state.grids.θl, state.grids.k)
    @. state.k = thermalConductivity(soil.tcparams, state.θw_k, state.θl_k, state.θm_k, state.θo_k)
    return nothing # ensure no allocation
end

function prognosticstep!(soil::Soil, heat::Heat{UT"J"}, state)
    δₓ = Δ(state.grids.T)
    # upper boundary
    state.dH[1] += let T₂=state.T[2], T₁=state.T[1], k=state.k[2], δ=δₓ[1];
        k*(T₂-T₁)/δ
    end
    # diffusion on non-boundary cells
    let T = state.T,
        k = (@view state.k[2:end-1])
        ∂H = (@view state.dH[2:end-1]);
        ∇²(T, δₓ, k, ∂H)
    end
    # lower boundary
    state.dH[end] += let T₂=state.T[end], T₁=state.T[end-1], k=state.k[end-1], δ=δₓ[end];
        -k*(T₂-T₁)/δ
    end
    return nothing # ensure no allocation
end

"""
Top interaction, constant temperature (Dirichlet) boundary condition.
"""
function interact!(top::Top, c::ConstantAirTemp, soil::Soil, heat::Heat{UT"J"}, stop, ssoil)
    δₓ = ssoil.grids.T[1] - ssoil.grids.k[1] # distance to surface
    ssoil.dH[1] += let Tair=c.value, Tsoil=ssoil.T[1], k=ssoil.k[1];
        -k*(Tsoil-Tair)/δₓ
    end
    return nothing # ensure no allocation
end

"""
Bottom interaction, constant geothermal heat flux (Neumann) boundary condition.
"""
function interact!(soil::Soil, heat::Heat{UT"J"}, bottom::Bottom, Qgeo::GeothermalHeatFlux, ssoil, sbot)
    ssoil.dH[end] += Qgeo.value
    return nothing # ensure no allocation
end
