abstract type RichardsEqFormulation end
struct SaturationBased <: RichardsEqFormulation end
struct PressureBased <: RichardsEqFormulation end
Base.@kwdef struct RichardsEq{Tform<:RichardsEqFormulation,Tswrc<:SWRCFunction,TΩ,Tsp} <: Hydrology.WaterFluxesImpl
    form::Tform = Pressure()
    swrc::Tswrc = VanGenuchten()
    Ω::TΩ = 1e-3 # smoothness for ice impedence factor
    advection_only::Bool = false
    sp::Tsp = nothing
end
"""
    impedencefactor(water::WaterFluxes{<:RichardsEq}, θw, θwi)

Impedence factor which represents the blockage of water-filled pores by ice (see Hansson et al. 2004 and Westermann et al. 2022).
"""
impedencefactor(water::WaterFluxes{<:RichardsEq}, θw, θwi) = 10^(-water.impl.Ω*(1 - θw/θwi))
swrc(water::WaterFluxes{<:RichardsEq}) = water.impl.swrc
Hydrology.saturation(soil::Soil, ::WaterFluxes, state, i) = porosity(soil, state, i)
@inline function Hydrology.waterice!(sub::SubSurface, water::WaterFluxes{<:RichardsEq{PressureBased}}, state)
    let swrc = swrc(water);
        @inbounds for i in 1:length(state.ψ)
            state.θsat[i] = saturation(sub, water, state, i)
            state.θwi[i], state.dθwidψ[i] = ∇(ψ -> swrc(ψ; θsat=state.θsat[i]), state.ψ[i])
            state.sat[i] = state.θwi[i] / state.θsat[i]
        end
    end
end
@inline function Hydrology.hydraulicconductivity!(sub::SubSurface, water::WaterFluxes{<:RichardsEq{Tform,<:VanGenuchten}}, state) where {Tform}
    kw_sat = kwsat(sub, water)
    vg = swrc(water)
    Δkw = Δ(state.grids.kw)
    @inbounds for i in 1:length(state.kwc)
        let θsat = saturation(sub, water, state, i),
            θw = state.θw[i],
            θwi = state.θwi[i],
            I_ice = impedencefactor(water, θw, θwi),
            n = vg.n;
            # van Genuchten formulation of hydraulic conductivity; see van Genuchten (1980) and Westermann et al. (2022).
            state.kwc[i] = kw_sat*I_ice*sqrt(θw/θsat)*(1 - (1 - (θw/θsat)^(n/(n+1)))^((n-1)/n))^2
        end
    end
    state.kw[1] = state.kwc[1]
    state.kw[end] = state.kwc[end]
    Numerics.harmonicmean!(@view(state.kw[2:end-1]), state.kwc, Δkw)
end
# CryoGrid methods
CryoGrid.variables(::RichardsEq{PressureBased}) = (
    Prognostic(:ψ, OnGrid(Cells), domain=0..1), # autmoatically generates dψ
    Diagnostic(:sat, OnGrid(Cells), domain=0..1), # saturation (diagnostic)
    Diagnostic(:θsat, OnGrid(Cells), domain=0..1), # maximum volumetric water content (saturation point)
    Diagnostic(:θwi, OnGrid(Cells), domain=0..1), # total water content (liquid + ice)
    Diagnostic(:θw, OnGrid(Cells), domain=0..1), # liquid water content
    Diagnostic(:dθwidψ, OnGrid(Cells), domain=0..Inf), # derivative of SWRC w.r.t matric potential
    Diagnostic(:dθwidt, OnGrid(Cells)), # time derivative of total water content
    Diagnostic(:jw, OnGrid(Edges), u"m/s"),
    Diagnostic(:kw, OnGrid(Edges), u"m/s", domain=0..Inf),
    Diagnostic(:kwc, OnGrid(Cells), u"m/s", domain=0..Inf),
)
function CryoGrid.initialcondition!(soil::Soil, water::WaterFluxes{<:RichardsEq{PressureBased}}, state)
    let invswrc = inv(swrc(water));
        # this will allocate but it's only `initialcondition!` so who cares
        @. state.θsat = saturation(soil, water, state, 1:length(state.θsat))
        # assumes that `sat` has been provided as an initial condition
        @. state.θwi = state.θsat*state.sat
        @. state.ψ = invswrc(state.θwi; θsat=state.θsat)
    end
end
function CryoGrid.diagnosticstep!(soil::Soil, ps::Coupled2{<:WaterFluxes,<:Heat}, state)
    water, heat = ps
    # reset fluxes
    Hydrology.resetfluxes!(soil, water, state)
    # first set water/ice diagnostics
    Hydrology.waterice!(soil, water, state)
    # then evaluate diagnosticstep! for Heat
    CryoGrid.diagnosticstep!(soil, heat, state)
    # then hydraulic conductivity
    Hydrology.hydraulicconductivity!(soil, water, state)
end
function CryoGrid.prognosticstep!(soil::Soil, water::WaterFluxes{<:RichardsEq{PressureBased}}, state)
    # downward advection due to gravity (the +1 in Richard's Equation)
    Hydrology.wateradvection!(soil, water, state)
    if water.impl.advection_only
        # compute divergernce using only advective fluxes
        Numerics.divergence!(state.dθwidt, state.jw, Δ(state.grids.kw))
    else
        # non-linear diffusion of water potential (fluxes are added to those from advection)
        Numerics.nonlineardiffusion!(state.dθwidt, state.jw, state.ψ, Δ(state.grids.ψ), state.kw, Δ(state.grids.kw))
    end
    # divide by specific moisture capacity (derivative of the water retention curve) dθwidψ
    @. state.dψ = state.dθwidt / state.dθwidψ
end
