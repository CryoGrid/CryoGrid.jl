module HeatConduction

import CryoGrid: SubSurfaceProcess, BoundaryStyle, Dirichlet, Neumann, BoundaryProcess, Layer, Top, Bottom, SubSurface, Callback
import CryoGrid: diagnosticstep!, prognosticstep!, interact!, initialcondition!, boundaryflux, boundaryvalue, variables, callbacks, criterion, affect!

using CryoGrid.Physics
using CryoGrid.Physics.Boundaries
using CryoGrid.Physics.Water: VanGenuchten
using CryoGrid.Numerics
using CryoGrid.Numerics: nonlineardiffusion!, harmonicmean!, harmonicmean, heaviside
using CryoGrid.Layers: Soil, porosity
using CryoGrid.Utils

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

export enthalpy, heatcapacity, heatcapacity!, thermalconductivity, thermalconductivity!
export heatconduction!
include("heat.jl")
export ConstantTemp, GeothermalHeatFlux, TemperatureGradient, NFactor, Damping
include("heat_bc.jl")
export SFCC, DallAmico, Westermann, McKenzie, SFCCNewtonSolver
include("soil/soilheat.jl")

end
