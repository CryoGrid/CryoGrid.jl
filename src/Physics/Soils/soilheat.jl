# === Thermal properties ===
# We use methods with optional index arguments `i` to allow for implementations both
# where these variables are treated as constants and as state variables.
# In the latter case, specializations should override only the index-free form
# and return a state vector instead of a scalar. The `getscalar` function will
# handle both the scalar and vector case!
# Functions for retrieving constituents and volumetric fractions
@inline HeatConduction.thermalconductivities(soil::Soil, heat::Heat) = (heat.prop.kw, heat.prop.ki, heat.prop.ka, soil.prop.km, soil.prop.ko)
@inline HeatConduction.heatcapacities(soil::Soil, heat::Heat) = (heat.prop.cw, heat.prop.ci, heat.prop.ca, soil.prop.cm, soil.prop.co)
@inline function Physics.volumetricfractions(soil::Soil, heat::Heat, state, i)
    return let θwi = waterice(soil, heat, state, i),
        θw = liquidwater(soil, heat, state, i),
        θm = mineral(soil, state, i),
        θo = organic(soil, state, i),
        θa = 1.0 - θwi - θm - θo,
        θi = θwi - θw;
        (θw, θi, θa, θm, θo)
    end
end

"""
Initial condition for heat conduction (all state configurations) on soil layer w/ SFCC.
"""
function CryoGrid.initialcondition!(soil::Soil, heat::Heat{<:SFCC,Enthalpy}, state)
    HeatConduction.initialcondition!(soil, heat, freezecurve(heat), state)
    CryoGrid.diagnosticstep!(soil, heat, state)
end
"""
Initial condition for heat conduction (all state configurations) on soil layer w/ free water freeze curve.
"""
function CryoGrid.initialcondition!(soil::Soil, heat::Heat{FreeWater}, state)
    L = heat.L
    # initialize liquid water content based on temperature
    @inbounds for i in 1:length(state.T)
        θwi = waterice(soil, heat, state, i)
        state.θw[i] = ifelse(state.T[i] > 0.0, θwi, 0.0)
        state.C[i] = heatcapacity(soil, heat, volumetricfractions(soil, heat, state, i)...)
        state.H[i] = enthalpy(state.T[i], state.C[i], L, state.θw[i])
    end
end
"""
    freezethaw!(soil::Soil, heat::Heat{<:SFCC,Temperature}, state)
    freezethaw!(soil::Soil, heat::Heat{<:SFCC,Enthalpy}, state)

Updates state variables according to the specified SFCC function and solver.
For heat conduction with enthalpy, evaluation of the inverse enthalpy function is performed using the given solver.
For heat conduction with temperature, we can simply evaluate the freeze curve to get C_eff, θw, and H.
"""
function HeatConduction.freezethaw!(soil::Soil, heat::Heat{<:SFCC,Temperature}, state)
    sfcc = freezecurve(heat)
    @inbounds @fastmath let L = heat.L,
        f = sfcc.f,
        f_args = sfccargs(f,soil,heat,state),
        f_res = ForwardDiff.DiffResult(zero(eltype(state.T)), zero(eltype(state.T)));
        for i in 1:length(state.T)
            T = state.T[i]
            f_argsᵢ = Utils.selectat(i, identity, f_args)
            θw, dθdT = ∇(T -> f(T, f_args...), T)
            state.θw[i] = θw
            state.dθdT[i] = dθdT
            state.C[i] = C = heatcapacity(soil, heat, volumetricfractions(soil, heat, state, i)...)
            state.dHdT[i] = C + (L + T*(heat.prop.cw - heat.prop.ci))*dθdT
            state.H[i] = enthalpy(T, C, L, θw)
        end
    end
end
function HeatConduction.freezethaw!(soil::Soil, heat::Heat{<:SFCC{F,SFCCNewtonSolver},Enthalpy}, state) where {F}
    sfcc = freezecurve(heat)
    f_args = sfccargs(sfcc.f, soil, heat, state)
    @inbounds for i in 1:length(state.H)
        let f_argsᵢ = Utils.selectat(i, identity, f_args),
            f = sfcc.f;
            # Evaluate inverse enthalpy function
            T = enthalpyinv(soil, heat, state, i)
            # Here we apply the recovered temperature to the state variables;
            # Since we perform iteration on untracked variables, we need to
            # recompute θw, C, and T here with the tracked variables.
            # Note that this results in one additional freeze curve function evaluation.
            state.θw[i] = f(T, f_argsᵢ...)
            dθdT = ForwardDiff.derivative(T -> f(T, f_argsᵢ...), T)
            let H = state.H[i],
                L = heat.L,
                cw = heat.prop.cw,
                ci = heat.prop.ci,
                θw = state.θw[i],
                (_, θi, θa, θm, θo) = volumetricfractions(soil, heat, state, i);
                state.C[i] = heatcapacity(soil, heat, θw, θi, θa, θm, θo)
                state.T[i] = (H - L*θw) / state.C[i]
                state.dHdT[i] = HeatConduction.C_eff(state.T[i], state.C[i], L, dθdT, cw, ci)
            end
        end
    end
end
function HeatConduction.freezethaw!(soil::Soil{<:HomogeneousCharacteristicFractions}, heat::Heat{<:SFCC{F,<:SFCCPreSolver},Enthalpy}, state) where {F}
    solver = freezecurve(heat).solver
    f = solver.cache.f
    ∇f = solver.cache.∇f
    @inbounds for i in 1:length(state.H)
        H = state.H[i]
        state.θw[i] = θw = f(H)
        state.C[i] = C = heatcapacity(soil, heat, volumetricfractions(soil, heat, state, i)...)
        state.T[i] = T = enthalpyinv(H, C, heat.L, θw)
        dθdTᵢ = ∇f(state.H[i])
        state.dHdT[i] = C + dθdTᵢ*(heat.L + T*(heat.prop.cw - heat.prop.ci))
    end
end
