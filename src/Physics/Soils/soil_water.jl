"""
    RREqForm

Base type for different formulations of Richardson-Richard's equation.
"""
abstract type RREqForm end
"""
    Saturation <: RREqForm

Sautration based form of Richards equation:

∇ₜθ = ∇ₓ[k(θ)[1 + ∇ₓΨ(θ)]]
"""
struct Saturation <: RREqForm end
"""
    Pressure <: RREqForm

Pressure head based form of Richards equation:

∂θ∂ψ⋅∇ₜΨ = ∇ₓ[k(θ)[1 + ∇ₓΨ]]
"""
struct Pressure <: RREqForm end

"""
    RichardsEq{Tform<:RichardsEqFormulation,Tswrc<:SWRC,Taux,TΩ} <: Hydrology.WaterFlow

The Richardson-Richards equation describes the flow of water in porous media under unsaturated condition.
"""
Base.@kwdef struct RichardsEq{Tform<:RREqForm,Tswrc<:SWRC,Taux,TΩ} <: Hydrology.WaterFlow
    form::Tform = Saturation()
    swrc::Tswrc = VanGenuchten()
    Ω::TΩ = 7 # scaling for ice impedence
    aux::Taux = nothing
end

swrc(water::WaterBalance{<:RichardsEq}) = water.flow.swrc

"""
    impedencefactor(water::WaterBalance{<:RichardsEq}, θw, θwi)

Impedence factor which represents the blockage of water-filled pores by ice (see Hansson et al. 2004 and Westermann et al. 2022).
"""
impedencefactor(water::WaterBalance{<:RichardsEq}, θw, θwi) = 10^(-water.flow.Ω*(1 - θw/θwi))

# Methods for Hydrology module

function Hydrology.hydraulicconductivity(soil::Soil, water::WaterBalance{<:RichardsEq{<:RREqForm,<:VanGenuchten}}, θw, θwi, θsat)
    let kw_sat = Hydrology.kwsat(soil, water),
        n = swrc(water).n,
        I_ice = impedencefactor(water, θw, θwi);
        # van Genuchten formulation of hydraulic conductivity; see van Genuchten (1980) and Westermann et al. (2022).
        # we use `complex` types here to permit illegal state values which may occur for adaptive solving schemes
        kw = abs(kw_sat*I_ice*sqrt(complex(θw/θsat))*(1 - complex(1 - complex(θw/θsat)^(n/(n+1)))^((n-1)/n))^2)
        return kw
    end
end

Hydrology.default_dtlim(::RichardsEq{Pressure}) = CryoGrid.MaxDelta(0.01u"m")
Hydrology.default_dtlim(::RichardsEq{Saturation}) = CryoGrid.MaxDelta(0.005)

Hydrology.maxwater(soil::Soil, ::WaterBalance, state, i) = porosity(soil, state, i)

function Hydrology.watercontent!(soil::Soil, water::WaterBalance{<:RichardsEq{Pressure}}, state)
    let swrc = swrc(water);
        @inbounds for i in 1:length(state.ψ₀)
            state.θsat[i] = Hydrology.maxwater(soil, water, state, i)
            state.ψ[i] = state.ψ₀[i] # initially set liquid pressure head to total water pressure head
            state.θwi[i], state.dθwidψ[i] = ∇(ψ -> swrc(ψ; θsat=state.θsat[i]), state.ψ[i])
            state.θw[i] = state.θwi[i] # initially set liquid water content to total water content (coupling with HeatBalance will overwrite this)
            state.sat[i] = state.θwi[i] / state.θsat[i]
        end
    end
end
function Hydrology.watercontent!(soil::Soil, water::WaterBalance{<:RichardsEq{Saturation}}, state)
    let f = swrc(water),
        f⁻¹ = inv(f),
        θres = f.vol.θres;
        @inbounds for i in 1:length(state.ψ₀)
            state.θsat[i] = θsat = Hydrology.maxwater(soil, water, state, i)
            state.θwi[i] = θwi = state.sat[i]*θsat
            state.θw[i] = state.θwi[i] # initially set liquid water content to total water content (coupling with HeatBalance will overwrite this)
            # this is a bit shady because we're allowing for incorrect/out-of-bounds values of θwi, but this is necessary
            # for solving schemes that might attempt to use illegal state values
            state.ψ₀[i] = f⁻¹(max(θres, min(θwi, θsat)); θsat)
            state.ψ[i] = state.ψ₀[i] # initially set liquid pressure head to total water pressure head
        end
    end
end

function Hydrology.waterprognostic!(::Soil, ::WaterBalance{<:RichardsEq{Saturation}}, state)
    @inbounds @. state.dsat = state.dθwi / state.θsat
    return nothing
end
function Hydrology.waterprognostic!(::Soil, ::WaterBalance{<:RichardsEq{Pressure}}, state)
    @inbounds @. state.dψ₀ = state.dθwi / state.∂θw∂ψ
    return nothing
end

function Hydrology.waterdiffusion!(::Soil, water::WaterBalance{<:RichardsEq}, state)
    # compute diffusive fluxes from pressure, if enabled
    Numerics.flux!(state.jw_v, state.ψ, Δ(cells(state.grid)), state.kw)
    return nothing
end

# CryoGrid methods
CryoGrid.variables(::RichardsEq{Pressure}) = (
    Prognostic(:ψ₀, OnGrid(Cells), domain=-Inf..0), # soil matric potential of water + ice
    Diagnostic(:ψ, OnGrid(Cells), domain=-Inf..0), # soil matric potential of unfrozen water
    Diagnostic(:sat, OnGrid(Cells), domain=0..1), # saturation (diagnostic)
    Diagnostic(:∂θw∂ψ, OnGrid(Cells), domain=0..Inf), # derivative of SWRC w.r.t matric potential
)

CryoGrid.variables(::RichardsEq{Saturation}) = (
    Prognostic(:sat, OnGrid(Cells), domain=0..1), # saturation
    Diagnostic(:ψ₀, OnGrid(Cells), domain=-Inf..0), # soil matric potential of water + ice
    Diagnostic(:ψ, OnGrid(Cells), domain=-Inf..0), # soil matric potential of unfrozen water
)

function CryoGrid.interact!(soil1::Soil, water1::WaterBalance{<:RichardsEq}, soil2::Soil, water2::WaterBalance{<:RichardsEq}, state1, state2)
    θw₁ = state1.θw[end]
    θw₂ = state2.θw[1]
    θwi₁ = state1.θwi[end]
    θwi₂ = state2.θwi[1]
    θsat₁ = state1.θsat[end]
    θsat₂ = state2.θsat[1]
    sat₁ = state1.sat[end]
    sat₂ = state2.sat[1]
    ψ₁ = state1.ψ[end]
    ψ₂ = state1.ψ[1]
    z₁ = CryoGrid.midpoint(soil1, state1, last)
    z₂ = CryoGrid.midpoint(soil2, state2, first)
    Δz₁ = CryoGrid.thickness(soil1, state1, last)
    Δz₂ = CryoGrid.thickness(soil2, state2, first)
    # take minimum water content from upper layer where water would drain from
    θmin₁ = Hydrology.minwater(soil1, water1, state1, lastindex(state1.θw))
    # hydraulic conductivity
    kwc₁ = state1.kwc[end]
    kwc₂ = state2.kwc[1]
    kw = state1.kw[end] = state2.kw[1] = harmonicmean(kwc₁, kwc₂, Δz₁, Δz₂)
    # flux over boundary = advective flux + diffusive pressure flux
    jw_adv = Hydrology.advectiveflux(θw₁, θmin₁, kw)
    jw_dif = -kw*(ψ₂ - ψ₁)/(z₂ - z₁)
    jw = (jw_adv + jw_dif)*state1.dt
    # setting both jw[end] on the upper layer and jw[1] on the lower layer is redundant since they refer to the same
    # element of the same underlying state array, but it's nice for clarity
    state1.jw[end] = state2.jw[1] = Hydrology.balanceflux(water1, water2, jw, θw₁, θw₂, θwi₁, θwi₂, θsat₁, θsat₂, sat₁, sat₂, Δz₁, Δz₂)
    return nothing
end

function CryoGrid.timestep(::Soil, water::WaterBalance{<:RichardsEq{Pressure},TET,<:CryoGrid.MaxDelta}, state) where {TET}
    dtmax = Inf
    @inbounds for i in 1:length(state.sat)
        dt = water.dtlim(state.dψ₀[i], state.ψ[i], state.t, -Inf, zero(eltype(state.ψ)))
        dt = isfinite(dt) ? dt : Inf # make sure it's +Inf
        dtmax = min(dtmax, dt)
    end
    return dtmax
end
