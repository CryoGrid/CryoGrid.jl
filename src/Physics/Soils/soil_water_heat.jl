default_sfccsolver(::Coupled(WaterBalance,HeatBalance)) = SFCCPreSolver(Solvers.SFCCPreSolverCacheND())
# Initialization
function CryoGrid.initialcondition!(soil::Soil, ps::Coupled(WaterBalance, HeatBalance), state)
    water, heat = ps
    CryoGrid.initialcondition!(soil, water, state)
    CryoGrid.initialcondition!(soil, heat, state)
end
function CryoGrid.initialcondition!(
    soil::Soil,
    ps::Coupled(WaterBalance{<:RichardsEq}, HeatBalance{<:SFCC}),
    state,
)
    water, heat = ps
    # initialize water
    CryoGrid.diagnosticstep!(soil, water, state)
    # initialize heat
    fc = heat.freezecurve
    solver = sfccsolver(soil)
    L = heat.prop.L
    hc = partial_heatcapacity(soil, heat)
    θsat = porosity(soil, state)
    @unpack ch_w, ch_i = thermalproperties(soil)
    FreezeCurves.Solvers.initialize!(solver, fc, hc; θsat)
    @inbounds for i in 1:length(state.T)
        fc_kwargsᵢ = sfcckwargs(fc, soil, heat, state, i)
        T = state.T[i]
        θw, ∂θw∂T = ∇(T -> fc(T, state.sat[i]; fc_kwargsᵢ...), T)
        state.θw[i] = θw
        state.C[i] = heatcapacity(soil, heat, volumetricfractions(soil, state, i)...)
        state.H[i] = enthalpy(state.T[i], state.C[i], L, state.θw[i])
        state.∂H∂T[i] = Heat.dHdT(T, state.C[i], L, ∂θw∂T, ch_w, ch_i)
    end
end
# Freeze/thaw dynamics
function Heat.freezethaw!(
    soil::Soil,
    ps::Coupled(WaterBalance{<:RichardsEq{TREqForm}}, HeatBalance{<:SFCC,THeatForm}),
    state
) where {TREqForm,THeatForm}
    water, heat = ps
    sfcc = heat.freezecurve
    swrc = FreezeCurves.swrc(sfcc.f)
    # helper function for computing temperature (inverse enthalpy, if necessary)
    _get_temperature(::Type{<:Temperature}, i) = state.T[i]
    _get_temperature(::Type{<:Enthalpy}, i) = enthalpyinv(soil, heat, state, i)
    L = heat.prop.L
    @unpack ch_w, ch_i = thermalproperties(soil)
    @inbounds @fastmath for i in 1:length(state.T)
        T = _get_temperature(THeatForm, i)
        ψ₀ = state.ψ₀[i]
        θtot = state.θwi[i]
        θsat = state.θsat[i]
        (ψ, ∂ψ∂T) = ∇(Tᵢ -> sfcc(Tᵢ, ψ₀, Val{:ψ}(); θtot, θsat), T)
        (θw, ∂θw∂ψ) = ∇(ψᵢ -> swrc(ψᵢ; θsat), ψ)
        ∂θw∂T = ∂θw∂ψ*∂ψ∂T
        C = Heat.heatcapacity(soil, heat, volumetricfractions(soil, state, i)...)
        ∂H∂T = Heat.dHdT(T, C, L, ∂θw∂T, ch_w, ch_i)
        state.∂θw∂T[i] = ∂θw∂T
        state.θw[i] = θw
        state.ψ[i] = ψ
        # compute dependent quantities
        state.C[i] = C
        state.∂H∂T[i] = ∂H∂T
        if TREqForm == Pressure
            state.∂θw∂ψ[i] = ∂θw∂ψ
        end
        if THeatForm == Temperature
            state.H[i] = Heat.enthalpy(T, C, L, θw)
        end
    end
    return nothing
end
function Heat.freezethaw!(
    soil::Soil,
    ps::Coupled(WaterBalance{<:BucketScheme}, HeatBalance{FreeWater}),
    state
)
    Heat.freezethaw!(soil, ps[2], state)
end