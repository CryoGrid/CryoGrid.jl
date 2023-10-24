# Freeze/thaw dynamics
Heat.freezethaw!(soil::Soil, ps::Coupled(WaterBalance, HeatBalance), state) = Heat.freezethaw!(soil, ps[2], state)
# special implementation of freezethaw! for the pressure-head  based form;
# this is necessary because of the need to compute ∂θw∂ψ
function Heat.freezethaw!(
    soil::Soil,
    ps::Coupled(WaterBalance{<:RichardsEq{Pressure}},HeatBalance{<:SFCC,THeatForm}),
    state
) where {THeatForm<:Heat.HeatOperator}
    water, heat = ps
    sfcc = heat.freezecurve
    swrc = FreezeCurves.swrc(sfcc)
    # helper function for computing temperature (inverse enthalpy, if necessary)
    _get_temperature(::Type{<:TemperatureBased}, i) = state.T[i]
    _get_temperature(::Type{<:EnthalpyBased}, i) = enthalpyinv(soil, heat, state, i)
    L = heat.prop.L
    @inbounds @fastmath for i in 1:length(state.T)
        @unpack ch_w, ch_i = thermalproperties(soil, state, i)
        T = state.T[i] = _get_temperature(THeatForm, i)
        sat = state.sat[i]
        θsat = state.θsat[i]
        ψ, ∂ψ∂T = ∇(Tᵢ -> sfcc(Tᵢ, sat, Val{:ψ}(); θsat), T)
        θw, ∂θw∂ψ = ∇(ψᵢ -> swrc(ψᵢ; θsat), ψ)
        ∂θw∂T = ∂θw∂ψ*∂ψ∂T
        state.θw[i] = θw
        state.ψ[i] = ψ
        C = Heat.heatcapacity(soil, heat, state, i)
        ∂H∂T = Heat.dHdT(T, C, L, ∂θw∂T, ch_w, ch_i)
        state.∂θw∂T[i] = ∂θw∂T
        # compute dependent quantities
        state.C[i] = C
        state.∂H∂T[i] = ∂H∂T
        state.∂θw∂ψ[i] = ∂θw∂ψ
        if THeatForm == TemperatureBased
            state.H[i] = Heat.enthalpy(T, C, L, θw)
        end
    end
    return nothing
end

# Initialization
function CryoGrid.initialcondition!(
    soil::Soil,
    ps::Coupled(WaterBalance{<:RichardsEq}, HeatBalance{<:SFCC}),
    state,
)
    water, heat = ps
    # initialize water
    CryoGrid.computediagnostic!(soil, water, state)
    # initialize heat
    fc = heat.freezecurve
    solver = initialize_sfccsolver!(soil, heat, state)
    @assert !isnothing(solver) "SFCC solver must be provided in HeatBalance operator. Check the model configuration."
    L = heat.prop.L
    @inbounds for i in 1:length(state.T)
        @unpack ch_w, ch_i = thermalproperties(soil, state, i)
        fc_kwargsᵢ = sfcckwargs(fc, soil, heat, state, i)
        T = state.T[i]
        θw, ∂θw∂T = ∇(T -> fc(T, state.sat[i]; fc_kwargsᵢ...), T)
        state.θw[i] = θw
        state.C[i] = heatcapacity(soil, heat, state, i)
        state.H[i] = enthalpy(state.T[i], state.C[i], L, state.θw[i])
        state.∂H∂T[i] = Heat.dHdT(T, state.C[i], L, ∂θw∂T, ch_w, ch_i)
    end
end

# Diagnostic step
function CryoGrid.computediagnostic!(sub::SubSurface, ps::Coupled(WaterBalance, HeatBalance), state)
    water, heat = ps
    # Compute water contents from current state
    Hydrology.watercontent!(sub, water, state)
    # HeatBalance diagnostics
    # Evaluate freeze/thaw processes
    Heat.freezethaw!(sub, state)
    # Update thermal conductivity
    Heat.thermalconductivity!(sub, heat, state)
    # then hydraulic conductivity (requires liquid water content from heat conduction)
    Hydrology.hydraulicconductivity!(sub, water, state)
end
