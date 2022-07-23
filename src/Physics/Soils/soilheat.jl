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
# Default implementations of variables for SFCCFunction
CryoGrid.variables(::Soil, ::Heat, f::SFCCFunction) = ()
CryoGrid.variables(::Soil, ::Heat, s::SFCCSolver) = ()
# Join the declared state variables of the SFCC function and the solver
CryoGrid.variables(soil::Soil, heat::Heat, sfcc::SFCC) = tuplejoin(CryoGrid.variables(soil, heat, sfcc.f), CryoGrid.variables(soil, heat, sfcc.solver))
"""
    sfcckwargs(f::SFCCFunction, soil::Soil, heat::Heat, state, i)

Builds a named tuple of values corresponding to each keyword arguments of the SFCCFunction `f`
which should be set according to the layer/process properties or state. The default implementation
sets only the total water content, θtot = θwi, and the saturated water content, θsat = θp.
"""
sfcckwargs(::SFCCFunction, soil::Soil, heat::Heat, state, i) = (
    θtot = waterice(soil, heat, state, i), # total water content    
    θsat = porosity(soil, state, i), # θ saturated = porosity
)
"""
Initial condition for heat conduction (all state configurations) on soil layer w/ SFCC.
"""
function CryoGrid.initialcondition!(soil::Soil, heat::Heat{<:SFCC}, state)
    fc = freezecurve(heat)
    L = heat.prop.L
    @inbounds for i in 1:length(state.T)
        fc_kwargsᵢ = sfcckwargs(fc.f, soil, heat, state, i)
        hc = partial(heatcapacity, liquidwater, soil, heat, state, i)
        if i == 1
            # TODO: this is currently only relevant for the pre-solver scheme and assumes that
            # the total water content is uniform throughout the layer and does not change over time.
            FreezeCurves.Solvers.initialize!(fc.solver, fc.f, hc; fc_kwargsᵢ...)
        end
        T = state.T[i]
        θw, dθdT = ∇(T -> fc(T; fc_kwargsᵢ...), T)
        state.θw[i] = θw
        state.C[i] = hc(θw)
        state.H[i] = enthalpy(state.T[i], state.C[i], L, state.θw[i])
        state.dHdT[i] = HeatConduction.C_eff(T, state.C[i], L, dθdT, heat.prop.cw, heat.prop.ci)
    end
end
"""
Initial condition for heat conduction (all state configurations) on soil layer w/ free water freeze curve.
"""
function CryoGrid.initialcondition!(soil::Soil, heat::Heat{FreeWater,Enthalpy}, state)
    L = heat.prop.L
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
    let L = heat.prop.L,
        f = sfcc.f;
        @inbounds @fastmath for i in 1:length(state.T)
            T = state.T[i]
            f_argsᵢ = sfcckwargs(f, soil, heat, state, i)
            θw, dθdT = ∇(T -> f(T; f_argsᵢ...), T)
            state.θw[i] = θw
            state.dθdT[i] = dθdT
            state.C[i] = C = heatcapacity(soil, heat, volumetricfractions(soil, heat, state, i)...)
            state.dHdT[i] = HeatConduction.C_eff(T, C, L, dθdT, heat.prop.cw, heat.prop.ci)
            state.H[i] = enthalpy(T, C, L, θw)
        end
    end
end
function HeatConduction.freezethaw!(soil::Soil, heat::Heat{<:SFCC,Enthalpy}, state)
    sfcc = freezecurve(heat)
    @inbounds for i in 1:length(state.H)
        @inbounds let H = state.H[i], # enthalpy
            L = heat.prop.L,
            cw = heat.prop.cw,
            ci = heat.prop.ci,
            θwi = waterice(soil, heat, state, i), # total water content
            T₀ = i > 1 ? state.T[i-1] : FreezeCurves.freewater(H, θwi, L), # initial guess for T
            hc = partial(heatcapacity, liquidwater, soil, heat, state, i),
            f = sfcc.f,
            f_kwargsᵢ = sfcckwargs(f, soil, heat, state, i),
            obj = FreezeCurves.SFCCInverseEnthalpyObjective(f, f_kwargsᵢ, hc, L, H);
            res = FreezeCurves.sfccsolve(obj, sfcc.solver, T₀, Val{true}())
            state.T[i] = res.T
            state.θw[i] = res.θw
            state.C[i] = res.C
            state.dHdT[i] = HeatConduction.C_eff(state.T[i], state.C[i], L, res.dθwdT, cw, ci)
        end
    end
end
function HeatConduction.enthalpyinv(soil::Soil, heat::Heat{<:SFCC,Enthalpy}, state, i)
    sfcc = freezecurve(heat)
    @inbounds let H = state.H[i], # enthalpy
        L = heat.prop.L, # latent heat of fusion of water
        θwi = waterice(soil, heat, state, i), # total water content
        T₀ = i > 1 ? state.T[i-1] : freewater(H, θwi, L),
        hc = partial(heatcapacity, liquidwater, soil, heat, state, i),
        f = sfcc.f,
        f_kwargsᵢ = sfcckwargs(f, soil, heat, state, i),
        obj = FreezeCurves.SFCCInverseEnthalpyObjective(f, f_kwargsᵢ, hc, L, H);
        T_sol = FreezeCurves.sfccsolve(obj, sfcc.solver, T₀, Val{false}())
        return T_sol
    end
end
