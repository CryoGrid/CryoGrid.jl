"""
    sfccheatcap(soil::Soil, heat::HeatBalance, state, i)

Constructs a function `heatcap(θw,θwi,θsat)` that comptues the heat capacity given
the current state for grid cell `i`.
"""
function sfccheatcap(soil::Soil, heat::HeatBalance, state, i)
    f = Heat.heat_capacity_function(heat.op)
    cs = Heat.heatcapacities(soil, state, i)
    function heatcap(θw, θwi, θsat)
        θi = θwi - θw
        θa = θsat - θwi
        θm = state.θm[i]
        θo = state.θo[i]
        return f(cs, (θw, θi, θa, θm, θo))
    end
end
function sfccheatcap(soil::Soil{<:MineralOrganic}, heat::HeatBalance, state, i)
    f = Heat.heat_capacity_function(heat.op)
    cs = Heat.heatcapacities(soil, state, i)
    function heatcap(θw, θwi, θsat)
        θi = θwi - θw
        θa = θsat - θwi
        θm = (1-soil.para.org)*(1-θsat)
        θo = soil.para.org*(1-θsat)
        return f(cs, (θw, θi, θa, θm, θo))
    end
end

"""
    initialize_sfccsolver!(::Soil, heat::HeatBalance, state)

Retrieves and initializes the SFCC solver for the given layer.
"""
initialize_sfccsolver!(soil::Soil, heat::HeatBalance, state) = fcsolver(soil)
function initialize_sfccsolver!(soil::Soil, heat::HeatBalance{<:SFCC,TOp}, state) where {TOp<:Heat.HeatOperator}
    fc = heat.freezecurve
    solver = fcsolver(soil)
    sat = saturation(soil, state)
    θsat = porosity(soil, state)
    hc = sfccheatcap(soil, heat, state, 1)
    if TOp <: EnthalpyBased
        @assert !isnothing(solver) "SFCC solver must be provided in HeatBalance operator. Check the model configuration."
        FreezeCurves.initialize!(solver, fc, hc; sat, θsat)
    end
    return solver
end

"""
    sfcckwargs(f::SFCC, soil::Soil, heat::HeatBalance, state, i)

Builds a named tuple of values corresponding to each keyword arguments of the SFCC `f`
which should be set according to the layer/process properties or state. The default implementation
sets only the saturated water content, θsat = porosity.
"""
sfcckwargs(::SFCC, soil::Soil, heat::HeatBalance, state, i) = (
    θsat = porosity(soil, state, i), # θ saturated = porosity
)

# Define volumetricfractions for Soil layer
function CryoGrid.volumetricfractions(soil::Soil, state, i)
    return let θwi = Hydrology.watercontent(soil, state, i),
        θw = state.θw[i],
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
function CryoGrid.initialcondition!(soil::Soil, heat::HeatBalance{<:SFCC}, state)
    L = heat.prop.L
    fc = heat.freezecurve
    initialize_sfccsolver!(soil, heat, state)
    @inbounds for i in 1:length(state.T)
        @unpack ch_w, ch_i = thermalproperties(soil, state, i)
        fc_kwargsᵢ = sfcckwargs(fc, soil, heat, state, i)
        sat = saturation(soil, state, i)
        T = state.T[i]
        θw, ∂θw∂T = ∇(T -> fc(T, sat; fc_kwargsᵢ...), T)
        state.θw[i] = θw
        state.C[i] = heatcapacity(soil, heat, state, i)
        state.H[i] = enthalpy(state.T[i], state.C[i], L, state.θw[i])
        state.∂H∂T[i] = Heat.dHdT(T, state.C[i], L, ∂θw∂T, ch_w, ch_i)
    end
end

"""
Initial condition for heat conduction (all state configurations) on soil layer w/ free water freeze curve.
"""
function CryoGrid.initialcondition!(soil::Soil, heat::HeatBalance{FreeWater}, state)
    L = heat.prop.L
    # initialize liquid water content based on temperature
    @inbounds for i in 1:length(state.T)
        θwi = Hydrology.watercontent(soil, state, i)
        state.θw[i] = ifelse(state.T[i] > 0.0, θwi, 0.0)
        state.C[i] = heatcapacity(soil, heat, state, i)
        state.H[i] = enthalpy(state.T[i], state.C[i], L, state.θw[i])
    end
end

"""
    freezethaw!(soil::Soil, heat::HeatBalance{<:SFCC,<:TemperatureBased}, state)
    freezethaw!(soil::Soil, heat::HeatBalance{<:SFCC,<:EnthalpyBased}, state)

Updates state variables according to the specified SFCC function and solver.
For heat conduction with enthalpy, evaluation of the inverse enthalpy function is performed using the given solver.
For heat conduction with temperature, we can simply evaluate the freeze curve to get dHdT, θw, and H.
"""
function Heat.freezethaw!(soil::Soil, heat::HeatBalance{<:SFCC,<:TemperatureBased}, state)
    sfcc = heat.freezecurve
    L = heat.prop.L
    @inbounds @fastmath for i in 1:length(state.T)
        T = state.T[i]
        @unpack ch_w, ch_i = Heat.thermalproperties(soil, state, i)
        f_argsᵢ = sfcckwargs(sfcc, soil, heat, state, i)
        θw, ∂θw∂T = ∇(T -> sfcc(T; f_argsᵢ...), T)
        state.θw[i] = θw
        state.∂θw∂T[i] = ∂θw∂T
        state.C[i] = C = heatcapacity(soil, heat, state, i)
        state.∂H∂T[i] = Heat.dHdT(T, C, L, ∂θw∂T, ch_w, ch_i)
        state.H[i] = enthalpy(T, C, L, θw)
    end
end

# freezethaw! implementation for enthalpy and implicit enthalpy formulations
function Heat.freezethaw!(soil::Soil, heat::HeatBalance{<:SFCC,<:EnthalpyBased}, state)
    sfcc = heat.freezecurve
    solver = fcsolver(soil)
    @inbounds for i in 1:length(state.H)
        let H = state.H[i], # enthalpy
            L = heat.prop.L,
            props = thermalproperties(soil, state, i),
            ch_w = props.ch_w,
            ch_i = props.ch_i,
            θwi = Hydrology.watercontent(soil, state, i),
            por = porosity(soil, state, i),
            sat = θwi / por,
            T₀ = i > 1 ? state.T[i-1] : H/ch_w, # initial guess for T
            f = sfcc,
            f_kwargsᵢ = sfcckwargs(f, soil, heat, state, i),
            hc = sfccheatcap(soil, heat, state, i),
            obj = FreezeCurves.SFCCInverseEnthalpyObjective(f, f_kwargsᵢ, hc, L, H, sat);
            res = FreezeCurves.sfccsolve(obj, solver, T₀, Val{true}())
            state.T[i] = res.T
            state.θw[i] = res.θw
            state.C[i] = res.C
            state.∂θw∂T[i] = res.∂θw∂T
            state.∂H∂T[i] = Heat.dHdT(state.T[i], state.C[i], L, res.∂θw∂T, ch_w, ch_i)
        end
    end
end

function Heat.enthalpyinv(soil::Soil, heat::HeatBalance{<:SFCC,<:EnthalpyBased}, state, i)
    sfcc = heat.freezecurve
    solver = fcsolver(soil)
    @inbounds let H = state.H[i], # enthalpy
        L = heat.prop.L, # latent heat of fusion of water
        props = thermalproperties(soil, state, i),
        ch_w = props.ch_w,
        θwi = Hydrology.watercontent(soil, state, i),
        por = porosity(soil, state, i),
        sat = θwi / por,
        T₀ = i > 1 ? state.T[i-1] : H/ch_w, # initial guess for T
        f = sfcc,
        f_kwargsᵢ = sfcckwargs(f, soil, heat, state, i),
        hc = sfccheatcap(soil, heat, state, i),
        obj = FreezeCurves.SFCCInverseEnthalpyObjective(f, f_kwargsᵢ, hc, L, H, sat);
        T_sol = FreezeCurves.sfccsolve(obj, solver, T₀, Val{false}())
        return T_sol
    end
end

# Freeze curve parameters;
# Since the freeze curve functions are specified in FreezeCurves.jl, we must (or rather should) provide
# CryoGrid.parameterize implementations for them here to specify parameter information.
CryoGrid.parameterize(f::PainterKarra) = PainterKarra(
    f.freezethaw,
    CryoGrid.parameterize(f.β, domain=OpenInterval(0,Inf), desc="Painter-Karra fitting parmaeter which controls the influence of temperature on the matric potential."),
    CryoGrid.parameterize(f.ω, domain=0..(1/f.β), desc="Painter-Karra fitting parameter which controls the depression of the melting point from saturation level."),
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
    CryoGrid.parameterize(f.saltconc, domain=Interval{:closed,:open}(0,Inf), desc="Assumed salt concentration when not defined as a state variable."), # salt concentration
    f.R,
    f.g,
    CryoGrid.parameterize(f.swrc),
)
CryoGrid.parameterize(f::McKenzie) = McKenzie(
    f.freezethaw,
    f.vol,
    CryoGrid.parameterize(f.γ, domain=OpenInterval(0,Inf)),
)
CryoGrid.parameterize(f::Westermann) = McKenzie(
    f.freezethaw,
    f.vol,
    CryoGrid.parameterize(f.δ, domain=OpenInterval(0,Inf)),
)
CryoGrid.parameterize(f::VanGenuchten) = VanGenuchten(
    f.vol,
    CryoGrid.parameterize(f.α, units=u"1/m", domain=OpenInterval(0,Inf), desc="van Genuchten α parameter which corresponds to the inverse of air entry suction."),
    CryoGrid.parameterize(f.n, domain=OpenInterval(1,Inf), desc="van Genuchten n parameter which controls the shape of the curve. Smaller values generally produce longer tailed curves."),
)
CryoGrid.parameterize(f::BrooksCorey) = BrooksCorey(
    f.vol,
    CryoGrid.parameterize(f.ψₛ, domain=OpenInterval(0,Inf), desc="Brooks-Corey suction parameter."),
    CryoGrid.parameterize(f.λ, domain=OpenInterval(0,Inf), desc="Brooks-Corey shape parameter."),
)
# do not parameterize default freeze curve properties or SFCC solvers
CryoGrid.parameterize(prop::FreezeCurves.SoilFreezeThawProperties) = prop
CryoGrid.parameterize(solver::FreezeCurves.SFCCSolver) = solver
