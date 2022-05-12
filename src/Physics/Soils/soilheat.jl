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
function liquidwater(::SubSurface, heat::Heat{<:SFCC,Temperature}, state, i)
    sfcc = freezecurve(heat)
    f_args = tuplejoin((state.T,), sfccargs(sfcc.f, soil, heat, state))
    f_argsᵢ = Utils.selectat(i, identity, f_args)
    return sfcc.f(f_argsᵢ...)
end
function liquidwater(sub::SubSurface, heat::Heat{<:SFCC,Enthalpy}, state, i)
    T = enthalpyinv(sub, heat, state, i)
    sfcc = freezecurve(heat)
    f_args = sfccargs(sfcc.f, soil, heat, state)
    f_argsᵢ = Utils.selectat(i, identity, f_args)
    return sfcc.f(T, f_argsᵢ...)
end
"""
    freezethaw!(soil::Soil, heat::Heat{<:SFCC,Temperature}, state)
    freezethaw!(soil::Soil, heat::Heat{<:SFCC,Enthalpy}, state)

Updates state variables according to the specified SFCC function and solver.
For heat conduction with enthalpy, evaluation of the inverse enthalpy function is performed using the given solver.
For heat conduction with temperature, we can simply evaluate the freeze curve to get C_eff, θl, and H.
"""
function freezethaw!(soil::Soil, heat::Heat{<:SFCC,Temperature}, state)
    sfcc = freezecurve(heat)
    @inbounds @fastmath let L = heat.L,
        f = sfcc.f,
        f_args = tuplejoin((state.T,), sfccargs(f,soil,heat,state));
        for i in 1:length(state.T)
            f_argsᵢ = Utils.selectat(i, identity, f_args)
            state.θl[i] = f(f_argsᵢ...)
            state.C[i] = heatcapacity(soil, heat, state, i)
            state.H[i] = enthalpy(state.T[i], state.C[i], L, state.θl[i])
            state.dHdT[i] = state.C[i] + (L + state.T[i]*(heat.prop.cw - heat.prop.ci))*∇f(f_argsᵢ)
        end
    end
end
function freezethaw!(soil::Soil, heat::Heat{<:SFCC{F,SFCCNewtonSolver},Enthalpy}, state) where {F}
    sfcc = freezecurve(heat)
    f_args = sfccargs(sfcc.f, soil, heat, state)
    @inbounds for i in 1:length(state.H)
        let f_argsᵢ = Utils.selectat(i, identity, f_args),
            f = sfcc.f;
            # Evaluate inverse enthalpy function
            T = enthalpyinv(soil, heat, state, i)
            # Here we apply the recovered temperature to the state variables;
            # Since we perform iteration on untracked variables, we need to
            # recompute θl, C, and T here with the tracked variables.
            # Note that this results in one additional freeze curve function evaluation.
            state.θl[i] = f(T, f_argsᵢ...)
            dθdT = ForwardDiff.derivative(T -> f(T, f_argsᵢ...), T)
            let θl = state.θl[i],
                H = state.H[i],
                L = heat.L,
                cw = heat.prop.cw,
                ci = heat.prop.ci,
                θw = totalwater(soil, state, i), # total water content
                θm = mineral(soil, state, i), # mineral content
                θo = organic(soil, state, i), # organic content
                θw = totalwater(soil, state, i);
                state.C[i] = heatcapacity(soil, heat, θw, θl, θm, θo)
                state.T[i] = (H - L*θl) / state.C[i]
                state.dHdT[i] = HeatConduction.C_eff(state.T[i], state.C[i], L, dθdT, cw, ci)
            end
        end
    end
end
function freezethaw!(soil::Soil{<:HomogeneousCharacteristicFractions}, heat::Heat{<:SFCC{F,<:SFCCPreSolver},Enthalpy}, state) where {F}
    solver = freezecurve(heat).solver
    f = solver.cache.f
    ∇f = solver.cache.∇f
    state.θl .= f.(state.H)
    heatcapacity!(soil, heat, state)
    @. state.T = (state.H - heat.L*state.θl) / state.C
    @inbounds for i in 1:length(state.H)
        dθdTᵢ = ∇f(state.H[i])
        state.dHdT[i] = state.C[i] + dθdTᵢ*(heat.L + state.T[i]*(heat.prop.cw - heat.prop.ci))
    end
end
