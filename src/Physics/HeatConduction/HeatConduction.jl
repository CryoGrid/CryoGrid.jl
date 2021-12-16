module HeatConduction

import CryoGrid: SubSurfaceProcess, BoundaryStyle, Dirichlet, Neumann, BoundaryProcess, Layer, Top, Bottom, SubSurface, Callback
import CryoGrid: diagnosticstep!, prognosticstep!, interact!, initialcondition!, boundaryflux, boundaryvalue, variables, callbacks, criterion, affect!
import CryoGrid.Layers: Soil, totalwater, porosity, mineral, organic

using CryoGrid.Physics
using CryoGrid.Physics.Boundaries
using CryoGrid.Physics.Water: VanGenuchten
using CryoGrid.Numerics
using CryoGrid.Numerics: nonlineardiffusion!, harmonicmean!, harmonicmean, heaviside
using CryoGrid.Utils

using Base: @propagate_inbounds
using DimensionalData
using IfElse
using Interpolations: Linear, Flat
using IntervalSets
using ModelParameters
using Parameters
using Unitful

import Flatten: @flattenable, flattenable

export Heat, TemperatureProfile
export FreeWater, FreezeCurve, freezecurve
export enthalpy, heatcapacity, heatcapacity!, thermalconductivity, thermalconductivity!

abstract type FreezeCurve end
struct FreeWater <: FreezeCurve end

TemperatureProfile(pairs::Pair{<:DistQuantity,<:TempQuantity}...) =
    Profile(map(p -> uconvert(u"m", p[1]) => uconvert(u"°C", p[2]),pairs)...)
    
abstract type HeatImpl end
struct Enthalpy <: HeatImpl end
struct Temperature <: HeatImpl end

@with_kw struct Heat{F<:FreezeCurve,S} <: SubSurfaceProcess
    ρ::Float"kg/m^3" = 1000.0xu"kg/m^3" #[kg/m^3]
    Lsl::Float"J/kg" = 334000.0xu"J/kg" #[J/kg] (latent heat of fusion)
    L::Float"J/m^3" = ρ*Lsl             #[J/m^3] (specific latent heat of fusion)
    freezecurve::F = FreeWater()        # freeze curve, defautls to free water fc
    sp::S = Enthalpy()                  # specialization
end
# convenience constructors for specifying prognostic variable as symbol
Heat(var::Union{Symbol,Tuple{Vararg{Symbol}}}; kwargs...) = Heat(Val{var}(); kwargs...)
Heat(::Val{:H}; kwargs...) = Heat(;sp=Enthalpy(), kwargs...)
Heat(::Val{:T}; kwargs...) = Heat(;sp=Temperature(), kwargs...)
freezecurve(heat::Heat) = heat.freezecurve
# Default implementation of `variables` for freeze curve
variables(::SubSurface, ::Heat, ::FreezeCurve) = ()
# Fallback (error) implementation for freeze curve
(fc::FreezeCurve)(sub::SubSurface, heat::Heat, state) = error("freeze curve $(typeof(fc)) not implemented for $(typeof(heat)) on layer $(typeof(sub))")
# Thermal properties (generic)
"""
    enthalpy(T, C, L, θ) = T*C + L*θ

Discrete enthalpy function on temperature, heat capacity, specific latent heat of fusion, and liquid water content.
"""
@inline enthalpy(T, C, L, θ) = T*C + L*θ
"""
    totalwater(sub::SubSurface, heat::Heat, state)
    totalwater(sub::SubSurface, heat::Heat, state, i)

Retrieves the total water content for the given layer at grid cell `i`, if specified.
Defaults to using the scalar total water content defined on layer `sub`.
"""
@inline totalwater(sub::SubSurface, heat::Heat, state) = totalwater(sub)
@inline totalwater(sub::SubSurface, heat::Heat, state, i) = Utils.getscalar(totalwater(sub, heat, state), i)
"""
    heatcapacity(sub::SubSurface, heat::Heat, state, i)

Computes the heat capacity at grid cell `i` for the given layer from the current state.
"""
heatcapacity(sub::SubSurface, heat::Heat, state, i) = error("heatcapacity not defined for $(typeof(heat)) on $(typeof(sub))")
"""
    heatcapacity!(sub::SubSurface, heat::Heat, state)

Computes the heat capacity for the given layer from the current state and stores the result in-place in the state variable `C`.
"""
@inline function heatcapacity!(sub::SubSurface, heat::Heat, state)
    @inbounds for i in 1:length(state.T)
        state.C[i] = heatcapacity(sub, heat, state, i)
    end
end
"""
    thermalconductivity(sub::SubSurface, heat::Heat, state, i)

Computes the thrmal conductivity at grid cell `i` for the given layer from the current state.
"""
thermalconductivity(sub::SubSurface, heat::Heat, state, i) = error("thermalconductivity not defined for $(typeof(heat)) on $(typeof(sub))")
"""
    thermalconductivity!(sub::SubSurface, heat::Heat, state)

Computes the thermal conductivity for the given layer from the current state and stores the result in-place in the state variable `C`.
"""
@inline function thermalconductivity!(sub::SubSurface, heat::Heat, state)
    @inbounds for i in 1:length(state.T)
        state.kc[i] = thermalconductivity(sub, heat, state, i)
    end
end

export heatconduction!
include("heat.jl")
export ConstantTemp, GeothermalHeatFlux, TemperatureGradient, NFactor, Damping
include("heat_bc.jl")
export SFCC, DallAmico, Westermann, McKenzie, SFCCNewtonSolver
include("soil/soilheat.jl")

end
