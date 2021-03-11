const TFreeWater = HeatParams{FreeWaterFC}

"""
    enthalpy(T,water,hc)

Enthalpy at temperature T with the given water content and heat capacity.
"""
function enthalpy(::Heat{UT"J",TFreeWater}, T, C, liquidWater, L)
    let θ = liquidWater; #[Vol. fraction]
        H = (T-273.15)*C + θ*L
    end
end

function enthalpyInv(::Heat{UT"J",TFreeWater}, H, C, totalWater, L)
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
function freezethaw(::Heat{UT"J",TFreeWater}, H, totalWater, L)
    let θ = max(1.0e-8, totalWater), #[Vol. fraction]
        Lθ = L*θ,
        I_t = H > Lθ,
        I_c = H > 0.0 && H <= Lθ;
        liquidfraction = I_c*(H/Lθ) + I_t
    end
end

export enthalpy, enthalpyInv, freezethaw

function diagnosticstep!(soil::Soil, heat::Heat{UT"J",TFreeWater}, state)
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

# Note for future use: harmonic mean of thermal conductivities
# for i=2:N-1
#     kn(i,1) = (dxp(i,1)/(2*dxn(i))*kp(i,1).^-1 + dxp(i-1,1)/(2*dxn(i))*kp(i-1).^-1).^-1;
#     ks(i,1) = (dxp(i,1)/(2*dxs(i))*kp(i,1).^-1 + dxp(i+1,1)/(2*dxs(i))*kp(i+1).^-1).^-1;
# end
