abstract type RichardsEqFormulation end
struct Saturation <: RichardsEqFormulation end
struct Pressure <: RichardsEqFormulation end
Base.@kwdef struct RichardsEq{Tform<:RichardsEqFormulation,Tswrc<:SWRCFunction,Tsp,Tkwsat,TΩ} <: Hydrology.WaterFlow
    form::Tform = Saturation()
    swrc::Tswrc = VanGenuchten()
    kw_sat::Tkwsat = Param(1e-5, domain=0..Inf, units=u"m/s")
    Ω::TΩ = 1e-3 # smoothness for ice impedence factor
    advection_only::Bool = false
    sp::Tsp = nothing
end
Hydrology.default_dtlim(::RichardsEq{Pressure}) = Physics.MaxDelta(0.01u"m")
Hydrology.default_dtlim(::RichardsEq{Saturation}) = Physics.MaxDelta(0.01)
Hydrology.kwsat(::Soil, water::WaterBalance{<:RichardsEq}) = water.flow.kw_sat
Hydrology.fieldcapacity(::Soil, water::WaterBalance{<:RichardsEq}) = 0.0
"""
    impedencefactor(water::WaterBalance{<:RichardsEq}, θw, θwi)

Impedence factor which represents the blockage of water-filled pores by ice (see Hansson et al. 2004 and Westermann et al. 2022).
"""
impedencefactor(water::WaterBalance{<:RichardsEq}, θw, θwi) = 10^(-water.flow.Ω*(1 - θw/θwi))
swrc(water::WaterBalance{<:RichardsEq}) = water.flow.swrc
Hydrology.saturation(soil::Soil, ::WaterBalance, state, i) = porosity(soil, state, i)
@inline function Hydrology.waterice!(sub::SubSurface, water::WaterBalance{<:RichardsEq{Pressure}}, state)
    let swrc = swrc(water);
        @inbounds for i in 1:length(state.ψ₀)
            state.θsat[i] = Hydrology.saturation(sub, water, state, i)
            state.ψ[i] = state.ψ₀[i] # initially set liquid pressure head to total water pressure head
            state.θwi[i], state.dθwidψ[i] = ∇(ψ -> swrc(ψ; θsat=state.θsat[i]), state.ψ[i])
            state.sat[i] = state.θwi[i] / state.θsat[i]
        end
    end
end
@inline function Hydrology.waterice!(sub::SubSurface, water::WaterBalance{<:RichardsEq{Saturation}}, state)
    let swrc⁻¹ = inv(swrc(water));
        @inbounds for i in 1:length(state.ψ₀)
            state.θsat[i] = Hydrology.saturation(sub, water, state, i)
            state.θwi[i] = state.sat[i]*state.θsat[i]
            state.ψ₀[i] = swrc⁻¹(state.θwi[i]; θsat=state.θsat[i])
            state.ψ[i] = state.ψ₀[i] # initially set liquid pressure head to total water pressure head
        end
    end
end
@inline function Hydrology.hydraulicconductivity!(sub::SubSurface, water::WaterBalance{<:RichardsEq{Tform,<:VanGenuchten}}, state) where {Tform}
    kw_sat = Hydrology.kwsat(sub, water)
    vg = swrc(water)
    Δkw = Δ(state.grids.kw)
    @inbounds for i in 1:length(state.kwc)
        let θsat = Hydrology.saturation(sub, water, state, i),
            θw = state.θw[i],
            θwi = state.θwi[i],
            I_ice = impedencefactor(water, θw, θwi),
            n = vg.n;
            # van Genuchten formulation of hydraulic conductivity; see van Genuchten (1980) and Westermann et al. (2022).
            state.kwc[i] = abs(kw_sat*I_ice*sqrt(θw/θsat)*(1 - (1 - (θw/θsat)^(n/(n+1)))^((n-1)/n))^2)
        end
    end
    state.kw[1] = state.kwc[1]
    state.kw[end] = state.kwc[end]
    Numerics.harmonicmean!(@view(state.kw[2:end-1]), state.kwc, Δkw)
end
function HeatConduction.freezethaw!(
    soil::Soil,
    ps::Coupled(WaterBalance{<:RichardsEq}, Heat{<:SFCC,Temperature}),
    state
)
    water, heat = ps
    sfcc = freezecurve(heat)
    let L = heat.prop.L,
        f = sfcc.f;
        @inbounds @fastmath for i in 1:length(state.T)
            T = state.T[i]
            ψ₀ = state.ψ₀[i]
            f_argsᵢ = sfcckwargs(f, soil, heat, state, i)
            # get derivative function (which will return results as ForwardDiff.Dual types)
            ∇f = ∇(T -> f(T, ψ₀; f_argsᵢ...))
            # unpack results
            θw, ψ, _ = ∇f(T, ψ₀)
            # extract derivaitve from water content
            dθwdT = ForwardDiff.partials(θw)[1]
            # exract values
            state.θw[i] = ForwardDiff.value(θw)
            state.ψ[i] = ForwardDiff.value(ψ)
            # compute dependent quantities
            state.dθwdT[i] = dθwdT
            state.C[i] = C = HeatConduction.heatcapacity(soil, heat, volumetricfractions(soil, heat, state, i)...)
            state.dHdT[i] = HeatConduction.C_eff(T, C, L, dθwdT, heat.prop.cw, heat.prop.ci)
            # derivaitve of enthalpy w.r.t water content when dθwdT > 0; otherwise set to zero to avoid Inf
            state.dHdθw[i] = IfElse.ifelse(dθwdT > 0, state.dHdT[i] / dθwdT, zero(dθwdT))
            # enthalpy
            state.H[i] = HeatConduction.enthalpy(T, C, L, θw)
        end
    end
    return nothing
end
# CryoGrid methods
CryoGrid.variables(soil::Soil, water::WaterBalance{<:RichardsEq{Pressure}}) = (
    Prognostic(:ψ₀, OnGrid(Cells), domain=-Inf..0), # soil matric potential of water + ice
    Diagnostic(:ψ, OnGrid(Cells), domain=-Inf..0), # soil matric potential of unfrozen water
    Diagnostic(:sat, OnGrid(Cells), domain=0..1), # saturation (diagnostic)
    Diagnostic(:dθwidψ, OnGrid(Cells), domain=0..Inf), # derivative of SWRC w.r.t matric potential
    CryoGrid.variables(soil)...,
    Hydrology.watervariables(water)...,
)
CryoGrid.variables(soil::Soil, water::WaterBalance{<:RichardsEq{Saturation}}) = (
    Prognostic(:sat, OnGrid(Cells), domain=0..1), # saturation
    Diagnostic(:ψ₀, OnGrid(Cells), domain=-Inf..0), # soil matric potential of water + ice
    Diagnostic(:ψ, OnGrid(Cells), domain=-Inf..0), # soil matric potential of unfrozen water
    CryoGrid.variables(soil)...,
    Hydrology.watervariables(water)...,
)
CryoGrid.variables(soil::Soil, ps::Coupled(WaterBalance, Heat)) = (
    Diagnostic(:dHdθw, OnGrid(Cells), domain=0..Inf),
    CryoGrid.variables(soil, ps[1])...,
)
function CryoGrid.prognosticstep!(soil::Soil, water::WaterBalance{<:RichardsEq{TForm}}, state) where {TForm}
    # downward advection due to gravity (the +1 in Richard's Equation)
    Hydrology.wateradvection!(soil, water, state)
    if !water.flow.advection_only
        # compute diffusive fluxes from pressure, if enabled
        Numerics.flux!(state.jw, state.ψ, Δ(state.grids.ψ), state.kw)
    end
    # compute divergence for water fluxes in all cells
    Numerics.divergence!(state.dθwidt, state.jw, Δ(state.grids.kw))
    # compute time derivative differently depending on Richard's equation formulation;
    # this if statement will be compiled away since TForm is known at compile time
    if TForm == Saturation
        # scale by max saturation (porosity) just like with bucket scheme
        @. state.dsat = state.dθwidt / state.θsat
    elseif TForm == Pressure
        # divide by specific moisture capacity (derivative of the water retention curve) dθwidψ
        @. state.dψ₀ = state.dθwidt / state.dθwidψ
    else
        error("$TForm not supported")
    end
    return nothing
end
function CryoGrid.timestep(::SubSurface, water::WaterBalance{<:RichardsEq{Pressure},<:Physics.MaxDelta}, state)
    dtmax = Inf
    @inbounds for i in 1:length(state.sat)
        dt = water.dtlim(state.dψ₀[i], state.ψ[i], state.t, -Inf, zero(eltype(state.ψ)))
        dt = isfinite(dt) ? dt : Inf # make sure it's +Inf
        dtmax = min(dtmax, dt)
    end
    return dtmax
end
function CryoGrid.timestep(::SubSurface, water::WaterBalance{<:RichardsEq{Saturation},<:Physics.MaxDelta}, state)
    dtmax = Inf
    @inbounds for i in 1:length(state.sat)
        # solve for dt in:
        # sat + dt*∂sat∂t = 1 if ∂sat∂t > 0
        # sat + dt*∂sat∂t = 0 if ∂sat∂t <= 0
        # sets the maximum timestep to the dt which would saturate the grid cell;
        # will be Inf when ∂sat∂t = 0
        dt = water.dtlim(state.dsat[i], state.sat[i], state.t, zero(state.t), one(state.t))
        dt = isfinite(dt) && dt > zero(dt) ? dt : Inf # make sure it's +Inf
        dtmax = min(dtmax, dt)
    end
    return dtmax
end
# Heat coupling
function CryoGrid.initialcondition!(soil::Soil, ps::Coupled2{<:WaterBalance,<:Heat}, state)
    CryoGrid.diagnosticstep!(soil, ps, state)
end
function CryoGrid.diagnosticstep!(soil::Soil, ps::Coupled2{<:WaterBalance,<:Heat}, state)
    water, heat = ps
    # reset fluxes
    Hydrology.resetfluxes!(soil, water, state)
    HeatConduction.resetfluxes!(sub, heat, state)
    # first set water/ice diagnostics
    Hydrology.waterice!(soil, water, state)
    Hydrology.liquidwater!(soil, water, state)
    # Evaluate freeze/thaw processes
    HeatConduction.freezethaw!(sub, ps, state)
    # Update thermal conductivity
    HeatConduction.thermalconductivity!(sub, heat, state)
    # then hydraulic conductivity
    Hydrology.hydraulicconductivity!(soil, water, state)
end