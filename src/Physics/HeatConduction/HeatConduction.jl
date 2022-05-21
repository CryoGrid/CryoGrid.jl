module HeatConduction

using CryoGrid: SubSurfaceProcess, BoundaryProcess, Dirichlet, Neumann, Layer, Top, Bottom, SubSurface
using CryoGrid.InputOutput: Forcing
using CryoGrid.Physics
using CryoGrid.Physics.Boundaries
using CryoGrid.Numerics
using CryoGrid.Numerics: nonlineardiffusion!, harmonicmean!, harmonicmean, heaviside
using CryoGrid.Utils

using Base: @propagate_inbounds, @kwdef
using IfElse
using ModelParameters
using Unitful

import CryoGrid
import CryoGrid: BoundaryStyle
import CryoGrid: diagnosticstep!, prognosticstep!, interact!, initialcondition!, boundaryflux, boundaryvalue, variables, thickness, midpoint
import CryoGrid.Physics

export Heat, TemperatureProfile
export FreeWater, FreezeCurve, freezecurve

abstract type FreezeCurve end
struct FreeWater <: FreezeCurve end

"""
    TemperatureProfile(pairs::Pair{<:Union{DistQuantity,Param},<:Union{TempQuantity,Param}}...)

Convenience constructor for `Numerics.Profile` which automatically converts temperature quantities.
"""
TemperatureProfile(pairs::Pair{<:Union{DistQuantity,Param},<:Union{TempQuantity,Param}}...) = Profile(map(p -> uconvert(u"m", p[1]) => uconvert(u"°C", p[2]), pairs))

"""
    HeatParameterization

Base type for different numerical formulations of two-phase heat diffusion.
"""
abstract type HeatParameterization end
struct Enthalpy <: HeatParameterization end
struct Temperature <: HeatParameterization end

ThermalProperties(
    consts=Physics.Constants();
    ρw = consts.ρw,
    Lf = consts.Lf,
    kw = Param(0.57, units=u"W/m/K"), # thermal conductivity of water [Hillel(1982)]
    ki = Param(2.2, units=u"W/m/K"), # thermal conductivity of ice [Hillel(1982)]
    ka = Param(0.025, units=u"W/m/K"), # air [Hillel(1982)]
    cw = Param(4.2e6, units=u"J/K/m^3"), # heat capacity of water
    ci = Param(1.9e6, units=u"J/K/m^3"), # heat capacity of ice
    ca = Param(0.00125e6, units=u"J/K/m^3"), # heat capacity of air
) = (; ρw, Lf, kw, ki, ka, cw, ci, ca)

@kwdef struct Heat{Tfc<:FreezeCurve,TPara<:HeatParameterization,Tinit,TProp,TL} <: SubSurfaceProcess
    para::TPara = Enthalpy()
    prop::TProp = ThermalProperties()
    L::TL = prop.ρw*prop.Lf # [J/m^3] (specific latent heat of fusion of water)
    freezecurve::Tfc = FreeWater() # freeze curve, defautls to free water fc
    init::Tinit = nothing # optional initialization scheme
end
# convenience constructors for specifying prognostic variable as symbol
Heat(var::Symbol; kwargs...) = Heat(Val{var}(); kwargs...)
Heat(::Val{:H}; kwargs...) = Heat(;para=Enthalpy(), kwargs...)
Heat(::Val{:T}; kwargs...) = Heat(;para=Temperature(), kwargs...)

# getter functions
thermalproperties(heat::Heat) = heat.prop
freezecurve(heat::Heat) = heat.freezecurve

# Default implementation of `variables` for freeze curve
variables(::SubSurface, ::Heat, ::FreezeCurve) = ()

export HeatBC, ConstantTemp, GeothermalHeatFlux, TemperatureGradient, NFactor
include("heat_bc.jl")

export heatconduction!, enthalpy, waterice, liquidwater, freezethaw!, heatcapacity, heatcapacity!, thermalconductivity, thermalconductivity!
include("heat.jl")

end
