abstract type RichardsEqFormulation end
struct SaturationBased <: RichardsEqFormulation end
struct PressureBased <: RichardsEqFormulation end
Base.@kwdef struct RichardsEq{Tform<:RichardsEqFormulation,Tswrc<:SWRCFunction,Tsp,Tkwsat,TΩ} <: Hydrology.WaterBalanceFormulation
    form::Tform = PressureBased()
    swrc::Tswrc = VanGenuchten()
    kw_sat::Tkwsat = Param(1e-5, domain=0..Inf, units=u"m/s")
    Ω::TΩ = 1e-3 # smoothness for ice impedence factor
    advection_only::Bool = false
    sp::Tsp = nothing
end
Hydrology.kwsat(::Soil, water::WaterBalance{<:RichardsEq}) = water.form.kw_sat
Hydrology.fieldcapacity(::Soil, water::WaterBalance{<:RichardsEq}) = 0.0
"""
    impedencefactor(water::WaterBalance{<:RichardsEq}, θw, θwi)

Impedence factor which represents the blockage of water-filled pores by ice (see Hansson et al. 2004 and Westermann et al. 2022).
"""
impedencefactor(water::WaterBalance{<:RichardsEq}, θw, θwi) = 10^(-water.form.Ω*(1 - θw/θwi))
swrc(water::WaterBalance{<:RichardsEq}) = water.form.swrc
Hydrology.saturation(soil::Soil, ::WaterBalance, state, i) = porosity(soil, state, i)
@inline function Hydrology.waterice!(sub::SubSurface, water::WaterBalance{<:RichardsEq{PressureBased}}, state)
    let swrc = swrc(water);
        @inbounds for i in 1:length(state.ψ)
            state.θsat[i] = Hydrology.saturation(sub, water, state, i)
            state.θwi[i], state.dθwidψ[i] = ∇(ψ -> swrc(ψ; θsat=state.θsat[i]), state.ψ[i])
            state.sat[i] = state.θwi[i] / state.θsat[i]
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
@inline function Hydrology.resetfluxes!(::SubSurface, water::WaterBalance{<:RichardsEq}, state)
    @. state.jw = zero(eltype(state.jw))
    @. state.dθwidt = zero(eltype(state.dθwidt))
end
# CryoGrid methods
CryoGrid.variables(::WaterBalance{<:RichardsEq{PressureBased}}) = (
    Prognostic(:ψ, OnGrid(Cells), domain=-Inf..0), # autmoatically generates dψ
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
function CryoGrid.initialcondition!(soil::Soil, water::WaterBalance{<:RichardsEq{PressureBased}}, state)
    let invswrc = inv(swrc(water));
        for i in 1:length(state.sat)
            state.θsat[i] = Hydrology.saturation(soil, water, state, i)
            # assumes sat has been provided as initial condition
            state.θwi[i] = state.θsat[i]*state.sat[i]
            state.ψ[i] = invswrc(state.θwi[i]; θsat=state.θsat[i])
            @assert isfinite(state.ψ[i])
        end
    end
end
function CryoGrid.diagnosticstep!(soil::Soil, ps::Coupled2{<:WaterBalance,<:Heat}, state)
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
function CryoGrid.prognosticstep!(soil::Soil, water::WaterBalance{<:RichardsEq{PressureBased}}, state)
    # downward advection due to gravity (the +1 in Richard's Equation)
    Hydrology.wateradvection!(soil, water, state)
    if water.form.advection_only
        # compute divergernce using only advective fluxes
        Numerics.divergence!(state.dθwidt, state.jw, Δ(state.grids.kw))
    else
        # non-linear diffusion of water potential (fluxes are added to those from advection)
        Numerics.nonlineardiffusion!(state.dθwidt, state.jw, state.ψ, Δ(state.grids.ψ), state.kw, Δ(state.grids.kw))
    end
    # divide by specific moisture capacity (derivative of the water retention curve) dθwidψ
    @. state.dψ = state.dθwidt / state.dθwidψ
end
# function CryoGrid.timestep(::SubSurface, water::WaterBalance{<:RichardsEq{PressureBased},<:Physics.DomainLimiter}, state)
#     dtmax = Inf
#     @inbounds for i in 1:length(state.sat)
#         dt = water.dtlim(state.dψ[i], state.ψ[i]; upper_limit=zero(eltype(state.sat)), upper_limit_factor=0.5)
#         dt = isfinite(dt) ? dt : Inf # make sure it's +Inf
#         dtmax = min(dtmax, dt)
#     end
#     return dtmax
# end