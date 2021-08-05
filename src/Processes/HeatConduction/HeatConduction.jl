module HeatConduction

using Unitful: Temperature
import CryoGrid.Interface: BoundaryStyle, diagnosticstep!, prognosticstep!, interact!, initialcondition!, variables

using ..Processes
using ..Processes.Boundaries
using ..Processes.Water: VanGenuchten
using CryoGrid.Forcings
using CryoGrid.Interface
using CryoGrid.Numerics
using CryoGrid.Numerics: nonlineardiffusion!
using CryoGrid.Layers: Soil, SoilParams
using CryoGrid.Utils

using DimensionalData
using IfElse
using Interpolations: Linear, Flat
using IntervalSets
using Parameters
using Unitful

export Heat, HeatParams, TempProfile
export FreeWater, FreezeCurve, freezecurve
export ConstantTemp, GeothermalHeatFlux, NFactor, TemperatureGradient
export SFCC, DallAmico, Westermann, McKenzie, SFCCNewtonSolver
export enthalpy, heatcapacity, heatcapacity!, thermalconductivity, thermalconductivity!
export heatconduction!, boundaryflux

abstract type FreezeCurve end
struct FreeWater <: FreezeCurve end

@with_kw struct HeatParams{F<:FreezeCurve,S} <: Params
    ρ::Float"kg/m^3" = 1000.0xu"kg/m^3" #[kg/m^3]
    Lsl::Float"J/kg" = 334000.0xu"J/kg" #[J/kg] (latent heat of fusion)
    L::Float"J/m^3" = ρ*Lsl             #[J/m^3] (specific latent heat of fusion)
    freezecurve::F = FreeWater()        # freeze curve, defautls to free water fc
    sp::S = nothing
end

"""
Alias and constructor for Profile specific to temperature.
"""
TempProfile(pairs::Pair{<:DistQuantity, <:TempQuantity}...) = Profile([d=>(uconvert(u"°C",T),) for (d,T) in pairs]...;names=(:T,))

struct Heat{U,F<:FreezeCurve,S} <: SubSurfaceProcess
    params::HeatParams{F,S}
    profile::Union{Nothing,<:DimArray{UFloat"°C"}}
    function Heat{var}(profile::TProfile=nothing; kwargs...) where {var,TProfile<:Union{Nothing,<:DimArray{UFloat"°C"}}}
        @assert var in [:H,(:Hₛ,:Hₗ)] "Invalid Heat prognostic variable: $var; must be one of :H, (:Hs,:Hl), or :T"
        params = HeatParams(;kwargs...)
        new{var,typeof(params.freezecurve),typeof(params.sp)}(params,profile)
    end
    function Heat{:T}(profile::TProfile=nothing; kwargs...) where {TProfile<:Union{Nothing,<:DimArray{UFloat"°C"}}}
        @assert :freezecurve in keys(kwargs) "Freeze curve must be specified for prognostic T heat configuration."
        @assert !(typeof(kwargs[:freezecurve]) <: FreeWater) "Free water freeze curve is not compatible with prognostic T."
        params = HeatParams(;kwargs...)
        new{:T,typeof(params.freezecurve),typeof(params.sp)}(params,profile)
    end
end

Base.show(io::IO, h::Heat{U,F,S}) where {U,F,S} = print(io, "Heat{$U,$F,$S}($(h.params))")

freezecurve(heat::Heat) = heat.params.freezecurve
enthalpy(T::Number"°C", C::Number"J/K/m^3", L::Number"J/m^3", θ::Real) = T*C + L*θ
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
variables(layer::SubSurface, heat::Heat{:H}) = (
    Prognostic(:H, Float"J/m^3", OnGrid(Cells)),
    Diagnostic(:T, Float"°C", OnGrid(Cells)),
    Diagnostic(:C, Float"J//K/m^3", OnGrid(Cells)),
    Diagnostic(:Ceff, Float"J/K/m^3", OnGrid(Cells)),
    Diagnostic(:k, Float"W/m/K", OnGrid(Edges)),
    Diagnostic(:kc, Float"W//m/K", OnGrid(Cells)),
    # add freeze curve variables (if any are present)
    variables(layer, heat, freezecurve(heat))...,
)
""" Variable definitions for heat conduction (partitioned enthalpy) on any subsurface layer. """
variables(layer::SubSurface, heat::Heat{(:Hₛ,:Hₗ)}) = (
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
    # add freeze curve variables (if any are present)
    variables(layer, heat, freezecurve(heat))...,
)
""" Variable definitions for heat conduction (temperature) on any subsurface layer. """
variables(layer::SubSurface, heat::Heat{:T}) = (
    Prognostic(:T, Float"°C", OnGrid(Cells)),
    Diagnostic(:H, Float"J/m^3", OnGrid(Cells)),
    Diagnostic(:dH, Float"J/s/m^3", OnGrid(Cells)),
    Diagnostic(:C, Float"J/K/m^3", OnGrid(Cells)),
    Diagnostic(:Ceff, Float"J/K/m^3", OnGrid(Cells)),
    Diagnostic(:k, Float"W/m/K", OnGrid(Edges)),
    Diagnostic(:kc, Float"W/m/K", OnGrid(Cells)),
    # add freeze curve variables (if any are present)
    variables(layer, heat, freezecurve(heat))...,
)
""" Initial condition for heat conduction (all state configurations) on any subsurface layer. """
function initialcondition!(layer::SubSurface, heat::Heat, state)
    interpolateprofile!(heat.profile, state)
    L = heat.params.L
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
    # Interpolate thermal conductivity to boundary grid
    regrid!(state.k, state.kc, state.grids.kc, state.grids.k, Linear(), Flat())
    # TODO: harmonic mean of thermal conductivities (in MATLAB code)
    # for i=2:N-1
    #     kn(i,1) = (dxp(i,1)/(2*dxn(i))*kp(i,1).^-1 + dxp(i-1,1)/(2*dxn(i))*kp(i-1).^-1).^-1;
    #     ks(i,1) = (dxp(i,1)/(2*dxs(i))*kp(i,1).^-1 + dxp(i+1,1)/(2*dxs(i))*kp(i+1).^-1).^-1;
    # end
    return nothing # ensure no allocation
end
""" Prognostic step for heat conduction (enthalpy) on subsurface layer. """
function prognosticstep!(::SubSurface, ::Heat{:H}, state)
    Δk = Δ(state.grids.k) # cell sizes
    ΔT = Δ(state.grids.T)
    # Diffusion on non-boundary cells
    heatconduction!(state.dH,state.T,ΔT,state.k,Δk)
end
""" Prognostic step for heat conduction (partitioned enthalpy) on subsurface layer."""
function prognosticstep!(::SubSurface, heat::Heat{(:Hₛ,:Hₗ)}, state)
    Δk = Δ(state.grids.k) # cell sizes
    ΔT = Δ(state.grids.T)
    # Diffusion on non-boundary cells
    heatconduction!(state.dH,state.T,ΔT,state.k,Δk)
    let L = heat.params.L;
        @. state.dHₛ = state.dH / (L/state.C*state.dθdT + 1)
        # This could also be expressed via a mass matrix with 1
        # in the upper right block diagonal. But this is easier.
        @. state.dHₗ = state.dH - state.dHₛ
    end
end
""" Prognostic step for heat conduction (temperature) on subsurface layer. """
function prognosticstep!(::SubSurface, ::Heat{:T}, state)
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
    @inbounds let a = Δ(ssub.grids.k)[1]
        bc(top,sub,h,stop,ssub)/a
    end
boundaryflux(::Neumann, bot::Bottom, bc::BoundaryProcess, sub::SubSurface, h::Heat, sbot, ssub) =
    @inbounds let a = Δ(ssub.grids.k)[end]
        bc(bot,sub,h,sbot,ssub)/a
    end
function boundaryflux(::Dirichlet, top::Top, bc::BoundaryProcess, sub::SubSurface, h::Heat, stop, ssub)
    Δk = Δ(ssub.grids.k)
    @inbounds let Tupper=bc(top,sub,h,stop,ssub),
        Tsub=ssub.T[1],
        k=ssub.k[1],
        a=Δk[1],
        δ=(Δk[1]/2); # distance to surface
        -k*(Tsub-Tupper)/δ/a
    end
end
function boundaryflux(::Dirichlet, bot::Bottom, bc::BoundaryProcess, sub::SubSurface, h::Heat, sbot, ssub)
    Δk = Δ(ssub.grids.k)
    @inbounds let Tlower=bc(bot,sub,h,sbot,ssub),
        Tsub=ssub.T[end],
        k=ssub.k[end],
        a=Δk[end],
        δ=(Δk[end]/2); # distance to surface
        -k*(Tsub-Tlower)/δ/a
    end
end
"""
Generic top interaction. Computes flux dH at top cell.
"""
function interact!(top::Top, bc::BoundaryProcess, sub::SubSurface, heat::Heat, stop, ssub)
    @inbounds ssub.dH[1] += boundaryflux(top, bc, sub, heat, stop, ssub)
    return nothing # ensure no allocation
end
"""
Generic bottom interaction. Computes flux dH at bottom cell.
"""
function interact!(sub::SubSurface, heat::Heat, bot::Bottom, bc::BoundaryProcess, ssub, sbot)
    @inbounds ssub.dH[end] += boundaryflux(bot, bc, sub, heat, sbot, ssub)
    return nothing # ensure no allocation
end
"""
Generic subsurface interaction. Computes flux dH at boundary between subsurface layers.
"""
function interact!(::SubSurface, ::Heat, ::SubSurface, ::Heat, s1, s2)
    # calculate heat flux between cells
    Qᵢ = @inbounds let k = (2*s1.k[end]*s2.k[1]) / (s1.k[end] + s2.k[1]), # harmonic mean of thermal conductivities
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
@inline function (fc::FreeWater)(layer::SubSurface, heat::Heat{:H}, state)
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
    L = heat.params.L
    @. state.θl = freezethaw(state.H, L, state.θw)*state.θw
    heatcapacity!(layer, heat, state) # update heat capacity, C
    @. state.T = enthalpyinv(state.H, state.C, L, state.θw)
end

# Default implementation of `variables` for freeze curve
variables(::SubSurface, ::Heat, ::FreezeCurve) = ()
# Fallback (error) implementation for freeze curve
(fc::FreezeCurve)(layer::SubSurface, heat::Heat, state) =
    error("freeze curve $(typeof(fc)) not implemented for $(typeof(heat)) on layer $(typeof(layer))")

include("soil/soilheat.jl")
include("heat_bc.jl")

end
