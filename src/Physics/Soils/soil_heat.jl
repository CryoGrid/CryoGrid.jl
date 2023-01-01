# === Thermal properties ===
# We use methods with optional index arguments `i` to allow for implementations both
# where these variables are treated as constants and as state variables.
# In the latter case, specializations should override only the index-free form
# and return a state vector instead of a scalar. The `getscalar` function will
# handle both the scalar and vector case!
@inline function Heat.thermalconductivities(soil::Soil)
    @unpack kh_w, kh_i, kh_a, kh_m, kh_o = thermalproperties(soil)
    return kh_w, kh_i, kh_a, kh_m, kh_o
end
@inline function Heat.heatcapacities(soil::Soil)
    @unpack hc_w, hc_i, hc_a, hc_m, hc_o = thermalproperties(soil)
    return hc_w, hc_i, hc_a, hc_m, hc_o
end
# Define volumetricfractions for Soil layer
@inline function Physics.volumetricfractions(soil::Soil, state, i)
    return let θwi = state.θwi[i],
        θw = state.θw[i],
        θm = mineral(soil, state, i),
        θo = organic(soil, state, i),
        θa = 1.0 - θwi - θm - θo,
        θi = θwi - θw;
        (θw, θi, θa, θm, θo)
    end
end
# Soil thermal properties
SoilThermalProperties(
    ::HomogeneousMixture;
    kh_o=0.25u"W/m/K", # organic [Hillel (1982)]
    kh_m=3.8u"W/m/K", # mineral [Hillel (1982)]
    hc_o=2.5e6u"J/K/m^3", # heat capacity organic
    hc_m=2.0e6u"J/K/m^3", # heat capacity mineral
    thermal_props...,
) = ThermalProperties(; kh_o, hc_o, kh_m, hc_m, thermal_props...)
# Soil properties for heat processes.
SoilProperties(para::HomogeneousMixture, ::HeatBalance; heat=SoilThermalProperties(para)) = SoilProperties(; heat)
"""
Gets the `ThermalProperties` for the given soil layer.
"""
Heat.thermalproperties(soil::Soil) = soil.prop.heat
"""
    sfcckwargs(f::SFCCFunction, soil::Soil, heat::HeatBalance, state, i)

Builds a named tuple of values corresponding to each keyword arguments of the SFCCFunction `f`
which should be set according to the layer/process properties or state. The default implementation
sets only the total water content, θtot = θwi, and the saturated water content, θsat = θp.
"""
sfcckwargs(::SFCCFunction, soil::Soil, heat::HeatBalance, state, i) = (
    θtot = state.θwi[i], # total water content    
    θsat = porosity(soil, state, i), # θ saturated = porosity
)
"""
Initial condition for heat conduction (all state configurations) on soil layer w/ SFCC.
"""
function CryoGrid.initialcondition!(soil::Soil, heat::HeatBalance{<:SFCC}, state)
    fc = freezecurve(heat)
    L = heat.prop.L
    @unpack hc_w, hc_i = thermalproperties(soil)
    @inbounds for i in 1:length(state.T)
        fc_kwargsᵢ = sfcckwargs(fc.f, soil, heat, state, i)
        hc = partial(heatcapacity, Val{:θw}(), soil, heat, state, i)
        if i == 1
            # TODO: this is currently only relevant for the pre-solver scheme and assumes that
            # the total water content is uniform throughout the layer and does not change over time.
            FreezeCurves.Solvers.initialize!(fc.solver, fc.f, hc; fc_kwargsᵢ...)
        end
        T = state.T[i]
        θw, ∂θw∂T = ∇(T -> fc(T; fc_kwargsᵢ...), T)
        state.θw[i] = θw
        state.C[i] = hc(θw)
        state.H[i] = enthalpy(state.T[i], state.C[i], L, state.θw[i])
        state.∂H∂T[i] = Heat.C_eff(T, state.C[i], L, ∂θw∂T, hc_w, hc_i)
    end
end
"""
Initial condition for heat conduction (all state configurations) on soil layer w/ free water freeze curve.
"""
function CryoGrid.initialcondition!(soil::Soil, heat::HeatBalance{FreeWater}, state)
    L = heat.prop.L
    # initialize liquid water content based on temperature
    @inbounds for i in 1:length(state.T)
        θwi = state.θwi[i]
        state.θw[i] = ifelse(state.T[i] > 0.0, θwi, 0.0)
        state.C[i] = heatcapacity(soil, heat, volumetricfractions(soil, state, i)...)
        state.H[i] = enthalpy(state.T[i], state.C[i], L, state.θw[i])
    end
end
"""
    freezethaw!(soil::Soil, heat::HeatBalance{<:SFCC,<:Temperature}, state)
    freezethaw!(soil::Soil, heat::HeatBalance{<:SFCC,<:Enthalpy}, state)

Updates state variables according to the specified SFCC function and solver.
For heat conduction with enthalpy, evaluation of the inverse enthalpy function is performed using the given solver.
For heat conduction with temperature, we can simply evaluate the freeze curve to get C_eff, θw, and H.
"""
function Heat.freezethaw!(soil::Soil, heat::HeatBalance{<:SFCC,<:Temperature}, state)
    sfcc = freezecurve(heat)
    f = sfcc.f
    L = heat.prop.L
    @unpack hc_w, hc_i = Heat.thermalproperties(soil)
    @inbounds @fastmath for i in 1:length(state.T)
        T = state.T[i]
        f_argsᵢ = sfcckwargs(f, soil, heat, state, i)
        θw, ∂θw∂T = ∇(T -> f(T; f_argsᵢ...), T)
        state.θw[i] = θw
        state.∂θw∂T[i] = ∂θw∂T
        state.C[i] = C = heatcapacity(soil, heat, volumetricfractions(soil, state, i)...)
        state.∂H∂T[i] = Heat.C_eff(T, C, L, ∂θw∂T, hc_w, hc_i)
        state.H[i] = enthalpy(T, C, L, θw)
    end
end
# freezethaw! implementation for enthalpy and implicit enthalpy formulations
function Heat.freezethaw!(soil::Soil, heat::HeatBalance{<:SFCC,<:Enthalpy}, state)
    sfcc = freezecurve(heat)
    @inbounds for i in 1:length(state.H)
        let H = state.H[i], # enthalpy
            L = heat.prop.L,
            props = thermalproperties(soil),
            hc_w = props.hc_w,
            hc_i = props.hc_i,
            θwi = state.θwi[i], # total water content
            T₀ = i > 1 ? state.T[i-1] : FreezeCurves.freewater(H, θwi, L), # initial guess for T
            hc = partial(heatcapacity, Val{:θw}(), soil, heat, state, i),
            f = sfcc.f,
            f_kwargsᵢ = sfcckwargs(f, soil, heat, state, i),
            obj = FreezeCurves.SFCCInverseEnthalpyObjective(f, f_kwargsᵢ, hc, L, H, nothing);
            res = FreezeCurves.sfccsolve(obj, sfcc.solver, T₀, Val{true}())
            state.T[i] = res.T
            state.θw[i] = res.θw
            state.C[i] = res.C
            state.∂H∂T[i] = Heat.C_eff(state.T[i], state.C[i], L, res.∂θw∂T, hc_w, hc_i)
        end
    end
end
function Heat.enthalpyinv(soil::Soil, heat::HeatBalance{<:SFCC,<:Enthalpy}, state, i)
    sfcc = freezecurve(heat)
    @inbounds let H = state.H[i], # enthalpy
        L = heat.prop.L, # latent heat of fusion of water
        θwi = state.θwi[i], # total water content
        hc = partial(Heat.heatcapacity, Val{:θw}(), soil, heat, state, i),
        T₀ = i > 1 ? state.T[i-1] : (H - L*θwi) / hc(θwi),
        f = sfcc.f,
        f_kwargsᵢ = sfcckwargs(f, soil, heat, state, i),
        obj = FreezeCurves.SFCCInverseEnthalpyObjective(f, f_kwargsᵢ, hc, L, H, nothing);
        T_sol = FreezeCurves.sfccsolve(obj, sfcc.solver, T₀, Val{false}())
        return T_sol
    end
end
# Freeze curve parameters;
# Since the freeze curve functions are specified in FreezeCurves.jl, we must (or rather should) provide
# CryoGrid.parameterize implementations for them here to specify parameter information.
CryoGrid.parameterize(sfcc::SFCC) = SFCC(
    CryoGrid.parameterize(sfcc.f),
    sfcc.solver,
)
CryoGrid.parameterize(f::PainterKarra) = PainterKarra(
    f.freezethaw,
    CryoGrid.parameterize(f.β, domain=OpenInterval(0,Inf), desc="Painter-Karra fitting parmaeter which controls the influence of temperature on the matric potential."),
    CryoGrid.parameterize(f.ω, domain=0..(1/β), desc="Painter-Karra fitting parameter which controls the depression of the melting point from saturation level."),
    f.g,
    CryoGrid.parameterize(f.swrc),
)
CryoGrid.parameterize(f::DallAmico) = DallAmico(
    f.freezethaw,
    f.g,
    CryoGrid.parameterize(f.swrc),
)
CryoGrid.parameterize(f::DallAmicoSalt) = DallAmicoSalt(
    f.freezethaw,
    CryoGrid.parameterize(f.saltconc, domain=Interval{:closed,:open}(0,Inf)), # salt concentration
    f.R,
    f.g,
    CryoGrid.parameterize(f.swrc),
)
CryoGrid.parameterize(f::McKenzie) = McKenzie(
    f.freezethaw,
    f.water,
    CryoGrid.parameterize(f.γ, domain=OpenInterval(0,Inf)),
)
CryoGrid.parameterize(f::Westermann) = McKenzie(
    f.freezethaw,
    f.water,
    CryoGrid.parameterize(f.δ, domain=OpenInterval(0,Inf)),
)
CryoGrid.parameterize(f::VanGenuchten) = VanGenuchten(
    f.water,
    CryoGrid.parameterize(f.α, domain=OpenInterval(0,Inf)),
    CryoGrid.parameterize(f.n, domain=OpenInterval(1,Inf)),
)
CryoGrid.parameterize(f::BrooksCorey) = BrooksCorey(
    f.water,
    CryoGrid.parameterize(f.ψₛ, domain=OpenInterval(0,Inf)),
    CryoGrid.parameterize(f.λ, domain=OpenInterval(0,Inf)),
)
# do not parameterize default freeze curve properties
CryoGrid.parameterize(prop::FreezeCurves.SoilFreezeThawProperties) = prop
CryoGrid.parameterize(prop::FreezeCurves.SoilWaterProperties) = prop
