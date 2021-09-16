module HeatConduction

import CryoGrid: SubSurfaceProcess, BoundaryStyle, Dirichlet, Neumann, BoundaryProcess, Layer, Top, Bottom, SubSurface, Boundary
import CryoGrid: diagnosticstep!, prognosticstep!, interact!, initialcondition!, variables

using ..Processes
using ..Processes.Boundaries
using ..Processes.Water: VanGenuchten
using CryoGrid.Numerics
using CryoGrid.Numerics: nonlineardiffusion!, harmonicmean!, harmonicmean, heaviside
using CryoGrid.Layers: Soil, θp, θw, θm, θo
using CryoGrid.Utils

using DimensionalData
using IfElse
using Interpolations: Linear, Flat
using IntervalSets
using ModelParameters
using Parameters
using Unitful

import Flatten: @flattenable, flattenable

export Heat, HeatParams, TemperatureProfile
export FreeWater, FreezeCurve, freezecurve
export ConstantTemp, GeothermalHeatFlux, NFactor, TemperatureGradient
export SFCC, DallAmico, Westermann, McKenzie, SFCCNewtonSolver
export enthalpy, heatcapacity, heatcapacity!, thermalconductivity, thermalconductivity!
export heatconduction!, boundaryflux

abstract type FreezeCurve end
struct FreeWater <: FreezeCurve end

TemperatureProfile(pairs::Pair{<:DistQuantity,<:TempQuantity}...) =
    Profile(map(p -> uconvert(u"m", p[1]) => uconvert(u"°C", p[2]),pairs)...)

DEFAULT_TEMP_PROFILE = TemperatureProfile(0.0u"m" => -1.0u"°C", 1000.0u"m" => -1.0u"°C") # uniform -1 degrees C

abstract type HeatVariable end
struct Enthalpy <: HeatVariable end
struct PartitionedEnthalpy <: HeatVariable end
struct Temperature <: HeatVariable end

@with_kw struct Heat{F<:FreezeCurve,S,N} <: SubSurfaceProcess
    initialT::Profile{N,typeof(1.0u"m"),typeof(1.0u"°C")} = DEFAULT_TEMP_PROFILE
    ρ::Float"kg/m^3" = 1000.0xu"kg/m^3" #[kg/m^3]
    Lsl::Float"J/kg" = 334000.0xu"J/kg" #[J/kg] (latent heat of fusion)
    L::Float"J/m^3" = ρ*Lsl             #[J/m^3] (specific latent heat of fusion)
    freezecurve::F = FreeWater()        # freeze curve, defautls to free water fc
    sp::S = Enthalpy()
end
# convenience constructors for specifying prognostic variable as symbol
Heat(var::Union{Symbol,Tuple{Vararg{Symbol}}}; kwargs...) = Heat(Val{var}(); kwargs...)
Heat(::Val{:H}; kwargs...) = Heat(;sp=Enthalpy(), kwargs...)
Heat(::Val{(:Hₛ,:Hₗ)}; kwargs...) = Heat(;sp=PartitionedEnthalpy(), kwargs...)
Heat(::Val{:T}; kwargs...) = Heat(;sp=Temperature(), kwargs...)

freezecurve(heat::Heat) = heat.freezecurve
enthalpy(T::Number"°C", C::Number"J/K/m^3", L::Number"J/m^3", θ::Real) = T*C + L*θ
totalwater(layer::SubSurface, heat::Heat, state) = state.θw
heatcapacity!(layer::SubSurface, heat::Heat, state) = error("heatcapacity not defined for $(typeof(heat)) on $(typeof(layer))")
thermalconductivity!(layer::SubSurface, heat::Heat, state) = error("thermalconductivity not defined for $(typeof(heat)) on $(typeof(layer))")
"""
    heatconduction!(∂H,T,ΔT,k,Δk)

1-D heat conduction/diffusion given T, k, and their deltas. Resulting enthalpy gradient is stored in ∂H.
Note that this function does not perform bounds checking. It is up to the user to ensure that all variables are
arrays of the correct length.
"""
function heatconduction!(∂H,T,ΔT,k,Δk)
    # upper boundary
    @inbounds ∂H[1] += let T₂=T[2],
        T₁=T[1],
        k=k[2],
        δ=ΔT[1],
        a=Δk[1];
        k*(T₂-T₁)/δ/a
    end
    # diffusion on non-boundary cells
    @inbounds let T = T,
        k = (@view k[2:end-1]),
        Δk = (@view Δk[2:end-1]),
        ∂H = (@view ∂H[2:end-1]);
        nonlineardiffusion!(∂H, T, ΔT, k, Δk)
    end
    # lower boundary
    @inbounds ∂H[end] += let T₂=T[end],
        T₁=T[end-1],
        k=k[end-1],
        δ=ΔT[end],
        a=Δk[end];
        -k*(T₂-T₁)/δ/a
    end
    return nothing
end
""" Variable definitions for heat conduction (enthalpy) on any subsurface layer. """
variables(layer::SubSurface, heat::Heat{<:FreezeCurve,Enthalpy}) = (
    Prognostic(:H, Float"J/m^3", OnGrid(Cells)),
    Diagnostic(:T, Float"°C", OnGrid(Cells)),
    Diagnostic(:C, Float"J//K/m^3", OnGrid(Cells)),
    Diagnostic(:Ceff, Float"J/K/m^3", OnGrid(Cells)),
    Diagnostic(:k, Float"W/m/K", OnGrid(Edges)),
    Diagnostic(:kc, Float"W//m/K", OnGrid(Cells)),
    Diagnostic(:θl, Float"1", OnGrid(Cells)),
    # add freeze curve variables (if any are present)
    variables(layer, heat, freezecurve(heat))...,
)
""" Variable definitions for heat conduction (partitioned enthalpy) on any subsurface layer. """
variables(layer::SubSurface, heat::Heat{<:FreezeCurve,PartitionedEnthalpy}) = (
    Prognostic(:Hₛ, Float"J/m^3", OnGrid(Cells)),
    Prognostic(:Hₗ, Float"J/m^3", OnGrid(Cells)),
    Diagnostic(:dH, Float"J/s/m^3", OnGrid(Cells)),
    Diagnostic(:H, Float"J", OnGrid(Cells)),
    Diagnostic(:T, Float"°C", OnGrid(Cells)),
    Diagnostic(:C, Float"J/K/m^3", OnGrid(Cells)),
    Diagnostic(:Ceff, Float"J/K/m^3", OnGrid(Cells)),
    Diagnostic(:dθdT, Float"m/m", OnGrid(Cells)),
    Diagnostic(:k, Float"W/m/K", OnGrid(Edges)),
    Diagnostic(:kc, Float"W/m/K", OnGrid(Cells)),
    Diagnostic(:θl, Float"1", OnGrid(Cells)),
    # add freeze curve variables (if any are present)
    variables(layer, heat, freezecurve(heat))...,
)
""" Variable definitions for heat conduction (temperature) on any subsurface layer. """
variables(layer::SubSurface, heat::Heat{<:FreezeCurve,Temperature}) = (
    Prognostic(:T, Float"°C", OnGrid(Cells)),
    Diagnostic(:H, Float"J/m^3", OnGrid(Cells)),
    Diagnostic(:dH, Float"J/s/m^3", OnGrid(Cells)),
    Diagnostic(:C, Float"J/K/m^3", OnGrid(Cells)),
    Diagnostic(:Ceff, Float"J/K/m^3", OnGrid(Cells)),
    Diagnostic(:k, Float"W/m/K", OnGrid(Edges)),
    Diagnostic(:kc, Float"W/m/K", OnGrid(Cells)),
    Diagnostic(:θl, Float"1", OnGrid(Cells)),
    # add freeze curve variables (if any are present)
    variables(layer, heat, freezecurve(heat))...,
)
""" Initial condition for heat conduction (all state configurations) on any subsurface layer. """
function initialcondition!(layer::SubSurface, heat::Heat, state)
    T₀ = profile2array(heat.initialT; names=(:T,))
    interpolateprofile!(T₀, state)
    L = heat.L
    heatcapacity!(layer, heat, state)
    @. state.H = enthalpy(state.T, state.C, L, state.θl)
end
""" Diagonstic step for heat conduction (all state configurations) on any subsurface layer. """
function diagnosticstep!(layer::SubSurface, heat::Heat, state)
    # Reset energy flux to zero; this is redundant when H is the prognostic variable
    # but necessary when it is not.
    @. state.dH = zero(eltype(state.dH))
    # Evaluate the freeze curve (updates T, C, and θl)
    fc! = freezecurve(heat);
    fc!(layer,heat,state)
    # Update thermal conductivity
    thermalconductivity!(layer, heat, state)
    # Harmonic mean of conductivities
    @inbounds let k = (@view state.k[2:end-1]),
        Δk = Δ(state.grids.k);
        harmonicmean!(k, state.kc, Δk)
    end
    return nothing # ensure no allocation
end
""" Prognostic step for heat conduction (enthalpy) on subsurface layer. """
function prognosticstep!(::SubSurface, ::Heat{<:FreezeCurve,Enthalpy}, state)
    Δk = Δ(state.grids.k) # cell sizes
    ΔT = Δ(state.grids.T)
    # Diffusion on non-boundary cells
    heatconduction!(state.dH,state.T,ΔT,state.k,Δk)
end
""" Prognostic step for heat conduction (partitioned enthalpy) on subsurface layer."""
function prognosticstep!(::SubSurface, heat::Heat{<:FreezeCurve,PartitionedEnthalpy}, state)
    Δk = Δ(state.grids.k) # cell sizes
    ΔT = Δ(state.grids.T)
    # Diffusion on non-boundary cells
    heatconduction!(state.dH,state.T,ΔT,state.k,Δk)
    let L = heat.L;
        @. state.dHₛ = state.dH / (L/state.C*state.dθdT + 1)
        # This could also be expressed via a mass matrix with 1
        # in the upper right block diagonal. But this is easier.
        @. state.dHₗ = state.dH - state.dHₛ
    end
end
""" Prognostic step for heat conduction (temperature) on subsurface layer. """
function prognosticstep!(::SubSurface, ::Heat{<:FreezeCurve,Temperature}, state)
    Δk = Δ(state.grids.k) # cell sizes
    ΔT = Δ(state.grids.T)
    # Diffusion on non-boundary cells
    heatconduction!(state.dH,state.T,ΔT,state.k,Δk)
    # Compute temperature flux by dividing by C_eff;
    # C_eff should be computed by the freeze curve.
    @inbounds @. state.dT = state.dH / state.Ceff
    return nothing
end
"""
    boundaryflux(boundary::Boundary, bc::BoundaryProcess, sub::SubSurface, h::Heat, sbound, ssub)

Computes the flux dH/dt at the given boundary. Calls boundaryflux(BoundaryStyle(B),...) to allow for generic
implementations by boundary condition type.
"""
boundaryflux(boundary::Boundary, bc::BoundaryProcess, sub::SubSurface, h::Heat, sbound, ssub) =
    boundaryflux(BoundaryStyle(bc), boundary, bc, sub, h, sbound, ssub)
boundaryflux(::Neumann, top::Top, bc::BoundaryProcess, sub::SubSurface, h::Heat, stop, ssub) =
    @inbounds let δ₀ = Δ(ssub.grids.k)[1]
        bc(top,sub,h,stop,ssub)/δ₀
    end
boundaryflux(::Neumann, bot::Bottom, bc::BoundaryProcess, sub::SubSurface, h::Heat, sbot, ssub) =
    @inbounds let δ₀ = Δ(ssub.grids.k)[end]
        bc(bot,sub,h,sbot,ssub)/δ₀
    end
function boundaryflux(::Dirichlet, top::Top, bc::BoundaryProcess, sub::SubSurface, h::Heat, stop, ssub)
    Δk = Δ(ssub.grids.k)
    @inbounds let Tupper=bc(top,sub,h,stop,ssub),
        Tsub=ssub.T[1],
        k=ssub.k[1],
        δ=Δk[1],
        δ₀=(Δk[1]/2); # distance to boundary
        -k*(Tsub-Tupper)/δ/δ₀
    end
end
function boundaryflux(::Dirichlet, bot::Bottom, bc::BoundaryProcess, sub::SubSurface, h::Heat, sbot, ssub)
    Δk = Δ(ssub.grids.k)
    @inbounds let Tlower=bc(bot,sub,h,sbot,ssub),
        Tsub=ssub.T[end],
        k=ssub.k[end],
        δ=Δk[1],
        δ₀=(Δk[1]/2); # distance to boundary
        -k*(Tsub-Tlower)/δ/δ₀
    end
end
"""
Generic top interaction. Computes flux dH at top cell.
"""
function interact!(top::Top, bc::BoundaryProcess, sub::SubSurface, heat::Heat, stop, ssub)
    # thermal conductivity at boundary
    # assumes (1) k has already been computed, (2) surface conductivity = cell conductivity
    @inbounds ssub.k[1] = ssub.k[2]
    # boundary flux
    @inbounds ssub.dH[1] += boundaryflux(top, bc, sub, heat, stop, ssub)
    return nothing # ensure no allocation
end
"""
Generic bottom interaction. Computes flux dH at bottom cell.
"""
function interact!(sub::SubSurface, heat::Heat, bot::Bottom, bc::BoundaryProcess, ssub, sbot)
    # thermal conductivity at boundary
    # assumes (1) k has already been computed, (2) bottom conductivity = cell conductivity
    @inbounds ssub.k[end] = ssub.k[end-1]
    # boundary flux
    @inbounds ssub.dH[end] += boundaryflux(bot, bc, sub, heat, sbot, ssub)
    return nothing # ensure no allocation
end
"""
Generic subsurface interaction. Computes flux dH at boundary between subsurface layers.
"""
function interact!(::SubSurface, ::Heat, ::SubSurface, ::Heat, s1, s2)
    # thermal conductivity between cells
    @inbounds let k₁ = s1.kc[end],
        k₂ = s2.kc[1],
        Δ₁ = Δ(s1.grids.k)[end],
        Δ₂ = Δ(s2.grids.k)[1];
        k = harmonicmean(k₁, k₂, Δ₁, Δ₂);
        s1.k[end] = s2.k[1] = k
    end
    # calculate heat flux between cells
    Qᵢ = @inbounds let k₁ = s1.k[end],
        k₂ = s2.k[1],
        k = k₁ = k₂,
        δ = s2.grids.T[1] - s1.grids.T[end];
        k*(s2.T[1] - s1.T[end]) / δ
    end
    # add fluxes scaled by grid cell size
    @inbounds s1.dH[end] += Qᵢ / Δ(s1.grids.k)[end]
    @inbounds s2.dH[1] += -Qᵢ / Δ(s2.grids.k)[1]
    return nothing # ensure no allocation
end
# Free water freeze curve
"""
Implementation of "free water" freeze curve for any subsurface layer. Assumes that
'state' contains at least temperature (T), enthalpy (H), heat capacity (C),
total water content (θw), and liquid water content (θl).
"""
@inline function (fc::FreeWater)(layer::SubSurface, heat::Heat{FreeWater,Enthalpy}, state)
    @inline function enthalpyinv(H, C, L, θtot)
        let θtot = max(1.0e-8,θtot),
            Lθ = L*θtot,
            I_t = H > Lθ,
            I_f = H <= 0.0;
            (I_t*(H-Lθ) + I_f*H)/C
        end
    end
    @inline function freezethaw(H, L, θtot)
        let θtot = max(1.0e-8,θtot),
            Lθ = L*θtot,
            I_t = H > Lθ,
            I_c = (H > 0.0) && (H <= Lθ);
            I_c*(H/Lθ) + I_t
        end
    end
    let L = heat.L,
        θw = totalwater(layer, heat, state);
        @. state.θl = freezethaw(state.H, L, θw)*θw
        heatcapacity!(layer, heat, state) # update heat capacity, C
        @. state.T = enthalpyinv(state.H, state.C, L, θw)
    end
    return nothing
end

# Default implementation of `variables` for freeze curve
variables(::SubSurface, ::Heat, ::FreezeCurve) = ()
# Fallback (error) implementation for freeze curve
(fc::FreezeCurve)(layer::SubSurface, heat::Heat, state) =
    error("freeze curve $(typeof(fc)) not implemented for $(typeof(heat)) on layer $(typeof(layer))")

include("soil/soilheat.jl")
include("heat_bc.jl")

end
