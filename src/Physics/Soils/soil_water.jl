# Soil hydrological properties
"""
Soil properites for water processes.
"""
SoilWaterProperties(
    ::HomogeneousMixture;
    kw_sat=1e-5u"m/s",
    hydraulic_props...,
) = HydraulicProperties(; kw_sat, hydraulic_props...)
SoilProperties(para::HomogeneousMixture, ::WaterBalance; water=SoilWaterProperties(para)) = SoilProperties(; water)
"""
Gets the `HydraulicProperties` for the given soil layer.
"""
Hydrology.hydraulicproperties(soil::Soil) = soil.prop.water
"""
Base type for different formulations of Richard's equation.
"""
abstract type RichardsEqFormulation end
struct Saturation <: RichardsEqFormulation end
struct Pressure <: RichardsEqFormulation end
Base.@kwdef struct RichardsEq{Tform<:RichardsEqFormulation,Tswrc<:SWRCFunction,Tsp,TΩ} <: Hydrology.WaterFlow
    form::Tform = Saturation()
    swrc::Tswrc = VanGenuchten()
    Ω::TΩ = 1e-3 # smoothness for ice impedence factor
    advection_only::Bool = false
    sp::Tsp = nothing
end
Hydrology.default_dtlim(::RichardsEq{Pressure}) = Physics.MaxDelta(0.01u"m")
Hydrology.default_dtlim(::RichardsEq{Saturation}) = Physics.MaxDelta(0.01)
Hydrology.fieldcapacity(::Soil, water::WaterBalance{<:RichardsEq}) = 0.0
"""
    impedencefactor(water::WaterBalance{<:RichardsEq}, θw, θwi)

Impedence factor which represents the blockage of water-filled pores by ice (see Hansson et al. 2004 and Westermann et al. 2022).
"""
impedencefactor(water::WaterBalance{<:RichardsEq}, θw, θwi) = 10^(-water.flow.Ω*(1 - θw/θwi))
swrc(water::WaterBalance{<:RichardsEq}) = water.flow.swrc
Hydrology.maxwater(soil::Soil, ::WaterBalance, state, i) = porosity(soil, state, i)
@inline function Hydrology.watercontent!(sub::SubSurface, water::WaterBalance{<:RichardsEq{Pressure}}, state)
    let swrc = swrc(water);
        @inbounds for i in 1:length(state.ψ₀)
            state.θsat[i] = Hydrology.maxwater(sub, water, state, i)
            state.ψ[i] = state.ψ₀[i] # initially set liquid pressure head to total water pressure head
            state.θwi[i], state.dθwidψ[i] = ∇(ψ -> swrc(ψ; θsat=state.θsat[i]), state.ψ[i])
            state.θw[i] = state.θwi[i] # initially set liquid water content to total water content (coupling with HeatBalance will overwrite this)
            state.sat[i] = state.θwi[i] / state.θsat[i]
        end
    end
end
@inline function Hydrology.watercontent!(sub::SubSurface, water::WaterBalance{<:RichardsEq{Saturation}}, state)
    let f = swrc(water),
        f⁻¹ = inv(f),
        θres = f.water.θres;
        @inbounds for i in 1:length(state.ψ₀)
            state.θsat[i] = θsat = Hydrology.maxwater(sub, water, state, i)
            state.θwi[i] = θwi = state.sat[i]*θsat
            # this is a bit shady because we're allowing for incorrect/out-of-bounds values of θwi, but this is necessary
            # for solving schemes that might attempt to use illegal state values
            state.ψ₀[i] = f⁻¹(max(θres, min(θwi, θsat)); θsat)
            state.ψ[i] = state.ψ₀[i] # initially set liquid pressure head to total water pressure head
        end
    end
end
@inline function Hydrology.hydraulicconductivity!(sub::SubSurface, water::WaterBalance{<:RichardsEq{Tform,<:VanGenuchten}}, state) where {Tform}
    kw_sat = Hydrology.kwsat(sub, water)
    vg = Soils.swrc(water)
    Δkw = Δ(state.grids.kw)
    @inbounds for i in 1:length(state.kwc)
        let θsat = Hydrology.maxwater(sub, water, state, i),
            θw = state.θw[i],
            θwi = state.θwi[i],
            I_ice = Soils.impedencefactor(water, θw, θwi),
            n = vg.n;
            # van Genuchten formulation of hydraulic conductivity; see van Genuchten (1980) and Westermann et al. (2022).
            # we use `complex` types here to permit illegal state values which may occur for adaptive solving schemes
            state.kwc[i] = abs(kw_sat*I_ice*sqrt(complex(θw/θsat))*(1 - complex(1 - complex(θw/θsat)^(n/(n+1)))^((n-1)/n))^2)
        end
    end
    state.kw[1] = state.kwc[1]
    state.kw[end] = state.kwc[end]
    Numerics.harmonicmean!(@view(state.kw[2:end-1]), state.kwc, Δkw)
end
function Hydrology.waterprognostic!(::Soil, ::WaterBalance{<:RichardsEq{Saturation}}, state)
    @inbounds @. state.∂sat∂t = state.∂θwi∂t / state.θsat
    return nothing
end
function Hydrology.waterprognostic!(::Soil, ::WaterBalance{<:RichardsEq{Pressure}}, state)
    @inbounds @. state.∂ψ₀∂t = state.∂θwi∂t / state.∂θw∂ψ
    return nothing
end
function Hydrology.waterdiffusion!(::Soil, water::WaterBalance{<:RichardsEq}, state)
    if !water.flow.advection_only
        # compute diffusive fluxes from pressure, if enabled
        Numerics.flux!(state.jw, state.ψ, Δ(state.grids.ψ), state.kw)
    end
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
function CryoGrid.initialcondition!(soil::Soil, water::WaterBalance, state)
    CryoGrid.diagnosticstep!(soil, water, state)
end
function CryoGrid.initialcondition!(soil::Soil, ps::Coupled(WaterBalance, HeatBalance), state)
    water, heat = ps
    CryoGrid.initialcondition!(soil, water, state)
    CryoGrid.initialcondition!(soil, heat, state)
end
function CryoGrid.interact!(sub1::SubSurface, water1::WaterBalance{<:RichardsEq}, sub2::SubSurface, water2::WaterBalance{<:RichardsEq}, state1, state2)
    θw₁ = state1.θw[end]
    ψ₁ = state1.ψ[end]
    ψ₂ = state2.ψ[1]
    θfc = Hydrology.fieldcapacity(sub1, water1) # take field capacity from upper layer where water would drain from
    kwc₁ = state1.kwc[end]
    kwc₂ = state2.kwc[1]
    δ₁ = CryoGrid.thickness(sub1, state1, last)
    δ₂ = CryoGrid.thickness(sub2, state2, first)
    z₁ = CryoGrid.midpoint(sub1, state1, last)
    z₂ = CryoGrid.midpoint(sub2, state2, first)
    kw = state1.kw[end] = state2.kw[1] = Numerics.harmonicmean(kwc₁, kwc₂, δ₁, δ₂)
    # flux over boundary = advective flux + diffusive pressure flux
    jw = Hydrology.advectiveflux(θw₁, θfc, kw) - kw*(ψ₂ - ψ₁)/(z₂ - z₁)
    # reduction factors
    r₁ = Hydrology.reductionfactor(water1, state1.sat[end])
    r₂ = Hydrology.reductionfactor(water2, state2.sat[1])
    # setting both jw[end] on the upper layer and jw[1] on the lower layer is redundant since they refer to the same
    # element of the same underlying state array, but it's nice for clarity
    state1.jw[end] = state2.jw[1] = jw*r₁*(jw < zero(jw)) + jw*r₂*(jw >= zero(jw))
    return nothing
end
function CryoGrid.timestep(::SubSurface, water::WaterBalance{<:RichardsEq{Pressure},TET,<:Physics.MaxDelta}, state) where {TET}
    dtmax = Inf
    @inbounds for i in 1:length(state.sat)
        dt = water.dtlim(state.∂ψ₀∂t[i], state.ψ[i], state.t, -Inf, zero(eltype(state.ψ)))
        dt = isfinite(dt) ? dt : Inf # make sure it's +Inf
        dtmax = min(dtmax, dt)
    end
    return dtmax
end
