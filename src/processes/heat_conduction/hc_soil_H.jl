"""
Soil heat conduction implementation for enthalpy (H) prognostic state.
"""

""" Defined variables for heat conduction (enthalpy) on soil layer. """
variables(soil::Soil, heat::Heat{UT"J"}) = (
    Prognostic(:H, Float"J", OnGrid(Cells)),
    Diagnostic(:T, Float"K", OnGrid(Cells)),
    Diagnostic(:C, Float"J/(K*m^3)", OnGrid(Cells)),
    Diagnostic(:k, Float"W/(m*K)", OnGrid(Edges)),
    Diagnostic(:kc, Float"W/(m*K)", OnGrid(Cells)),
)

# convenience constants for HeatParams parametric types
const FreeWaterFC = HeatParams{FreeWater}

"""
    enthalpy(T,water,hc)

Enthalpy at temperature T with the given water content and heat capacity.
"""
function enthalpy(::Heat{UT"J",FreeWaterFC}, T, C, liquidWater, L)
    let θ = liquidWater; #[Vol. fraction]
        H = (T-273.15)*C + θ*L
    end
end

function enthalpyInv(::Heat{UT"J",FreeWaterFC}, H, C, totalWater, L)
    let θ = max(1.0e-8, totalWater), #[Vol. fraction]
        Lθ = L*θ,
        # indicator variables for thawed and frozen states respectively
        I_t = H > Lθ,
        I_f = H <= 0.0;
        T = (I_t*(H-Lθ) + I_f*H)/C + 273.15
    end
end

"""
Phase change with linear freeze curve. Assumes diagnostic liquid water variable. Should *not* be used with prognostic
water variable.
"""
function freezethaw(::Heat{UT"J",FreeWaterFC}, H, totalWater, L)
    let θ = max(1.0e-8, totalWater), #[Vol. fraction]
        Lθ = L*θ,
        I_t = H > Lθ,
        I_c = H > 0.0 && H <= Lθ;
        liquidfraction = I_c*(H/Lθ) + I_t
    end
end

export enthalpy, enthalpyInv, freezethaw

function initialcondition!(soil::Soil, heat::Heat{UT"J"}, state)
    interpolateprofile!(heat.profile, state)
    L = ρ(heat)*Lsl(heat)
    @. state.C = heatCapacity(soil.hcparams, state.θw, state.θl, state.θm, state.θo)
    @. state.H = enthalpy(heat, state.T, state.C, state.θl, L)
    # k lies on the boundary grid, so we have to regrid the soil properties
    @. state.kc = thermalConductivity(soil.tcparams, state.θw, state.θl, state.θm, state.θo)
    regrid!(state.k, state.kc, state.grids.kc, state.grids.k, Linear(), Flat())
end

function diagnosticstep!(soil::Soil, heat::Heat{UT"J",FreeWaterFC}, state)
    let ρ = heat.params.ρ,
        Lsl = heat.params.Lsl,
        L = ρ*Lsl; #[J/m^3];
        @. state.T = enthalpyInv(heat, state.H, state.C, state.θw, L)
        @. state.θl = freezethaw(heat, state.H, state.θw, L) # overwrite θl temporarily to avoid allocation
        @. state.θl = state.θw*state.θl
    end
    @. state.C = heatCapacity(soil.hcparams, state.θw, state.θl, state.θm, state.θo)
    @. state.kc = thermalConductivity(soil.tcparams, state.θw, state.θl, state.θm, state.θo)
    regrid!(state.k, state.kc, state.grids.kc, state.grids.k, Linear(), Flat())
    return nothing # ensure no allocation
end

function prognosticstep!(soil::Soil, heat::Heat{UT"J"}, state)
    Δk = Δ(state.grids.k) # cell sizes
    ΔT = Δ(state.grids.T)
    # upper boundary
    state.dH[1] += let T₂=state.T[2],
        T₁=state.T[1],
        k=state.k[2],
        δ=ΔT[1],
        m=Δk[1];
        k*(T₂-T₁)/δ/m
    end
    # diffusion on non-boundary cells
    let T = state.T,
        k = (@view state.k[2:end-1]),
        ∂H = (@view state.dH[2:end-1]);
        ∇²(T, ΔT, k, ∂H)
    end
    # lower boundary
    state.dH[end] += let T₂=state.T[end],
        T₁=state.T[end-1],
        k=state.k[end-1],
        δ=ΔT[end],
        m=Δk[end];
        -k*(T₂-T₁)/δ/m
    end
    return nothing # ensure no allocation
end

"""
Top interaction, constant temperature (Dirichlet) boundary condition.
"""
function interact!(top::Top, c::ConstantAirTemp, soil::Soil, heat::Heat{UT"J"}, stop, ssoil)
    Δk = Δ(ssoil.grids.k)
    ssoil.dH[1] += let Tair=c.value,
        Tsoil=ssoil.T[1],
        k=ssoil.k[1],
        m=Δk[1],
        δ=(Δk[1]/2); # distance to surface
        -k*(Tsoil-Tair)/δ/m
    end
    return nothing # ensure no allocation
end


"""
Bottom interaction, constant geothermal heat flux (Neumann) boundary condition.
"""
function interact!(soil::Soil, heat::Heat{UT"J"}, bottom::Bottom, Qgeo::GeothermalHeatFlux, ssoil, sbot)
    Δk = Δ(ssoil.grids.k)
    ssoil.dH[end] += Qgeo.value/Δk[end]
    return nothing # ensure no allocation
end

# Note for future use: harmonic mean of thermal conductivities
# for i=2:N-1
#     kn(i,1) = (dxp(i,1)/(2*dxn(i))*kp(i,1).^-1 + dxp(i-1,1)/(2*dxn(i))*kp(i-1).^-1).^-1;
#     ks(i,1) = (dxp(i,1)/(2*dxs(i))*kp(i,1).^-1 + dxp(i+1,1)/(2*dxs(i))*kp(i+1).^-1).^-1;
# end
