"""
    sfccheatcap(soil::Soil, state, i)

Constructs a function `heatcap(θw,θwi,θsat)` that comptues the heat capacity given
the current state for grid cell `i`.
"""
function sfccheatcap(soil::Soil, state, i)
    f = thermalproperties(soil).heatcap
    cs = Heat.heatcapacities(soil, state, i)
    function heatcap(θw, θwi, θsat)
        θi = θwi - θw
        θa = θsat - θwi
        θm = mineral(soil, state, i)*(1-θsat)
        θo = organic(soil, state, i)*(1-θsat)
        return f(cs, (θw, θi, θa, θm, θo))
    end
end

"""
    initialize_sfccsolver!([::FreezeCurve], ::Soil, state)

Retrieves and initializes the SFCC solver for the given layer.
"""
initialize_sfccsolver!(soil::Soil, state) = initialize_sfccsolver!(freezecurve(soil), soil, state)
initialize_sfccsolver!(::FreezeCurve, soil::Soil, state) = fcsolver(soil)
function initialize_sfccsolver!(sfcc::SFCC, soil::Soil, state)
    solver = fcsolver(soil)
    sat = saturation(soil, state)
    θsat = porosity(soil, state)
    hc = sfccheatcap(soil, state, 1)
    if isa(soil.heat, HeatBalance{<:EnthalpyBased})
        @assert !isnothing(solver) "SFCC solver must be provided in HeatBalance operator. Check the model configuration."
        FreezeCurves.initialize!(solver, sfcc, hc; sat, θsat)
    end
    return solver
end

function initialize_soil_heat!(::FreeWater, soil::Soil, state)
    L = soil.heat.prop.L
    # initialize liquid water content based on temperature
    @inbounds for i in 1:length(state.T)
        θwi = Hydrology.watercontent(soil, state, i)
        state.θw[i] = ifelse(state.T[i] > 0.0, θwi, 0.0)
        state.C[i] = heatcapacity(soil, state, i)
        state.H[i] = enthalpy(state.T[i], state.C[i], L, state.θw[i])
    end
end

function initialize_soil_heat!(sfcc::SFCC, soil::Soil, state)
    L = soil.heat.prop.L
    initialize_sfccsolver!(soil, state)
    @inbounds for i in 1:length(state.T)
        @unpack ch_w, ch_i = thermalproperties(soil, state, i)
        fc_kwargsᵢ = sfcckwargs(sfcc, soil, state, i)
        sat = saturation(soil, state, i)
        T = state.T[i]
        θw, ∂θw∂T = ∇(T -> sfcc(T, sat; fc_kwargsᵢ...), T)
        state.θw[i] = θw
        state.C[i] = heatcapacity(soil, state, i)
        state.H[i] = enthalpy(state.T[i], state.C[i], L, state.θw[i])
        state.∂H∂T[i] = Heat.dHdT(T, state.C[i], L, ∂θw∂T, ch_w, ch_i)
    end
end

"""
    sfcckwargs(f::SFCC, soil::Soil, state, i)

Builds a named tuple of values corresponding to each keyword arguments of the SFCC `f`
which should be set according to the layer/process properties or state. The default implementation
sets only the saturated water content, θsat = porosity.
"""
sfcckwargs(::SFCC, soil::Soil, state, i) = (
    θsat = porosity(soil, state, i), # θ saturated = porosity
)

function soil_volumetricfractions(soil::Soil, state, i)
    return let θwi = Hydrology.watercontent(soil, state, i),
        θw = state.θw[i],
        θm = mineral(soil, state, i),
        θo = organic(soil, state, i),
        θa = 1.0 - θwi - θm - θo,
        θi = θwi - θw;
        (θw, θi, θa, θm, θo)
    end
end

# Define volumetricfractions for Soil layer
CryoGrid.volumetricfractions(soil::Soil, state, i) = soil_volumetricfractions(soil, state, i)

CryoGrid.initialcondition!(soil::Soil, heat::HeatBalance, state) = initialize_soil_heat!(freezecurve(soil), soil, state)

"""
    freezethaw!(sfcc::SFCC, soil::Soil, heat::HeatBalance{<:TemperatureBased}, state)
    freezethaw!(sfcc::SFCC, soil::Soil, heat::HeatBalance{<:EnthalpyBased}, state)

Updates state variables according to the specified SFCC function and solver.
For heat conduction with enthalpy, evaluation of the inverse enthalpy function is performed using the given solver.
For heat conduction with temperature, we can simply evaluate the freeze curve to get dHdT, θw, and H.
"""
function Heat.freezethaw!(sfcc::SFCC, soil::Soil, heat::HeatBalance{<:TemperatureBased}, state)
    L = heat.prop.L
    @inbounds @fastmath for i in 1:length(state.T)
        T = state.T[i]
        @unpack ch_w, ch_i = Heat.thermalproperties(soil, state, i)
        f_argsᵢ = sfcckwargs(sfcc, soil, state, i)
        θw, ∂θw∂T = ∇(T -> sfcc(T; f_argsᵢ...), T)
        state.θw[i] = θw
        state.∂θw∂T[i] = ∂θw∂T
        state.C[i] = C = heatcapacity(soil, state, i)
        state.∂H∂T[i] = Heat.dHdT(T, C, L, ∂θw∂T, ch_w, ch_i)
        state.H[i] = enthalpy(T, C, L, θw)
    end
end

# freezethaw! implementation for enthalpy and implicit enthalpy formulations
function Heat.freezethaw!(sfcc::SFCC, soil::Soil, heat::HeatBalance{<:EnthalpyBased}, state)
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
            f_kwargsᵢ = sfcckwargs(f, soil, state, i),
            hc = sfccheatcap(soil, state, i),
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

function Heat.enthalpyinv(sfcc::SFCC, soil::Soil, heat::HeatBalance{<:EnthalpyBased}, state, i)
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
        f_kwargsᵢ = sfcckwargs(f, soil, state, i),
        hc = sfccheatcap(soil, state, i),
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
