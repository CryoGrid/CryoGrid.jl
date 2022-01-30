module HeatConduction

import CryoGrid: SubSurfaceProcess, BoundaryStyle, Dirichlet, Neumann, BoundaryProcess, Layer, Top, Bottom, SubSurface, Callback
import CryoGrid: diagnosticstep!, prognosticstep!, interact!, initialcondition!, boundaryflux, boundaryvalue, variables, callbacks, criterion, affect!

using CryoGrid.Physics
using CryoGrid.Physics.Boundaries
using CryoGrid.Numerics
using CryoGrid.Numerics: nonlineardiffusion!, harmonicmean!, harmonicmean, heaviside
using CryoGrid.Utils

using Base: @propagate_inbounds
using IfElse
using ModelParameters
using Parameters
using SimulationLogs
using Unitful

export Heat, TemperatureProfile
export FreeWater, FreezeCurve, freezecurve

abstract type FreezeCurve end
struct FreeWater <: FreezeCurve end

"""
    TemperatureProfile(pairs::Pair{<:DistQuantity,<:Union{TempQuantity,Param}}...)

Convenience constructor for `Numerics.Profile` which automatically converts temperature quantities.
"""
TemperatureProfile(pairs::Pair{<:DistQuantity,<:Union{TempQuantity,Param}}...) = Profile(map(p -> uconvert(u"m", p[1]) => uconvert(u"°C", p[2]), pairs)...)

"""
    HeatImpl

Base type for different numerical formulations of two-phase heat diffusion.
"""
abstract type HeatImpl end
struct Enthalpy <: HeatImpl end
struct Temperature <: HeatImpl end

@with_kw struct Heat{Tfc<:FreezeCurve,Tsp} <: SubSurfaceProcess
    ρ::Float"kg/m^3" = 1000.0xu"kg/m^3" #[kg/m^3]
    Lsl::Float"J/kg" = 334000.0xu"J/kg" #[J/kg] (latent heat of fusion)
    L::Float"J/m^3" = ρ*Lsl #[J/m^3] (specific latent heat of fusion)
    freezecurve::Tfc = FreeWater() # freeze curve, defautls to free water fc
    sp::Tsp = Enthalpy() # specialization
end
# convenience constructors for specifying prognostic variable as symbol
Heat(var::Symbol; kwargs...) = Heat(Val{var}(); kwargs...)
Heat(::Val{:H}; kwargs...) = Heat(;sp=Enthalpy(), kwargs...)
Heat(::Val{:T}; kwargs...) = Heat(;sp=Temperature(), kwargs...)
freezecurve(heat::Heat) = heat.freezecurve
# Default implementation of `variables` for freeze curve
variables(::SubSurface, ::Heat, ::FreezeCurve) = ()
# Fallback (error) implementation for freeze curve
(fc::FreezeCurve)(sub::SubSurface, heat::Heat, state) = error("freeze curve $(typeof(fc)) not implemented for $(typeof(heat)) on layer $(typeof(sub))")

export heatconduction!, enthalpy, totalwater, liquidwater, heatcapacity, heatcapacity!, thermalconductivity, thermalconductivity!
include("heat.jl")
export ConstantTemp, GeothermalHeatFlux, TemperatureGradient, NFactor, Damping
include("heat_bc.jl")

end
