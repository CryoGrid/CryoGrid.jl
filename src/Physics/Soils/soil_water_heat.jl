# Initialization
# Water/heat coupling
function CryoGrid.initialcondition!(sub::SubSurface, ps::Coupled(WaterBalance, HeatBalance), state)
    CryoGrid.diagnosticstep!(sub, ps, state)
end
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
    solver = Heat.fcsolver(heat)
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

# Diagnostic step
function CryoGrid.diagnosticstep!(sub::SubSurface, ps::Coupled(WaterBalance, HeatBalance), state)
    water, heat = ps
    # Reset fluxes
    Hydrology.resetfluxes!(sub, water, state)
    # Compute water contents from current state
    Hydrology.watercontent!(sub, water, state)
    # HeatBalance diagnostics
    Heat.diagnosticstep!(sub, heat, state)
    # then hydraulic conductivity (requires liquid water content from heat conduction)
    Hydrology.hydraulicconductivity!(sub, water, state)
end

# Prognostic step
function CryoGrid.prognosticstep!(sub::SubSurface, ps::Coupled(WaterBalance, HeatBalance), state)
    water, heat = ps
    CryoGrid.prognosticstep!(sub, water, state)
    L = heat.prop.L
    @unpack ch_w, ch_i = thermalproperties(sub)
    # heat flux due to change in water content
    # @. state.∂H∂t += state.∂θwi∂t*(state.T*(ch_w - ch_i) + L)
    CryoGrid.prognosticstep!(sub, heat, state)
end

# Freeze/thaw dynamics
Heat.freezethaw!(soil::Soil, ps::Coupled(WaterBalance, HeatBalance), state) = Heat.freezethaw!(soil, ps[2], state)
function Heat.freezethaw!(
    soil::Soil,
    ps::Coupled2{<:WaterBalance{<:RichardsEq{TREqForm}},<:HeatBalance{<:SFCC,THeatForm}},
    state
) where {TREqForm,THeatForm<:Heat.HeatOperator}
    water, heat = ps
    sfcc = heat.freezecurve
    swrc = FreezeCurves.swrc(sfcc)
    # helper function for computing temperature (inverse enthalpy, if necessary)
    _get_temperature(::Type{<:Temperature}, i) = state.T[i]
    _get_temperature(::Type{<:Enthalpy}, i) = enthalpyinv(soil, heat, state, i)
    L = heat.prop.L
    @unpack ch_w, ch_i = thermalproperties(soil)
    @inbounds @fastmath for i in 1:length(state.T)
        T = state.T[i] = _get_temperature(THeatForm, i)
        sat = state.sat[i]
        θsat = state.θsat[i]
        ψ, ∂ψ∂T = ∇(Tᵢ -> sfcc(Tᵢ, sat, Val{:ψ}(); θsat), T)
        θw, ∂θw∂ψ = ∇(ψᵢ -> swrc(ψᵢ; θsat), ψ)
        ∂θw∂T = ∂θw∂ψ*∂ψ∂T
        state.θw[i] = θw
        state.ψ[i] = ψ
        C = Heat.heatcapacity(soil, heat, volumetricfractions(soil, state, i)...)
        ∂H∂T = Heat.dHdT(T, C, L, ∂θw∂T, ch_w, ch_i)
        state.∂θw∂T[i] = ∂θw∂T
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
