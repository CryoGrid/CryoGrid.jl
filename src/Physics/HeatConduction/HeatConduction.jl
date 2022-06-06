module HeatConduction

using CryoGrid
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

abstract type StepLimiter end
@kwdef struct CFL <: StepLimiter
    fallback_dt::Float64 = 60.0 # fallback dt [s]
end

ThermalProperties(
    consts=Physics.Constants();
    ρw = consts.ρw,
    Lf = consts.Lf,
    L = consts.ρw*consts.Lf,
    kw = Param(0.57, units=u"W/m/K"), # thermal conductivity of water [Hillel(1982)]
    ki = Param(2.2, units=u"W/m/K"), # thermal conductivity of ice [Hillel(1982)]
    ka = Param(0.025, units=u"W/m/K"), # air [Hillel(1982)]
    cw = Param(4.2e6, units=u"J/K/m^3"), # heat capacity of water
    ci = Param(1.9e6, units=u"J/K/m^3"), # heat capacity of ice
    ca = Param(0.00125e6, units=u"J/K/m^3"), # heat capacity of air
) = (; ρw, Lf, L, kw, ki, ka, cw, ci, ca)

struct Heat{Tfc<:FreezeCurve,TPara<:HeatParameterization,Tdt,Tinit,TProp} <: SubSurfaceProcess
    para::TPara
    prop::TProp
    freezecurve::Tfc
    dtlim::Tdt  # timestep limiter
    init::Tinit # optional initialization scheme
end
# convenience constructors for specifying prognostic variable as symbol
Heat(var::Symbol=:H; kwargs...) = Heat(Val{var}(); kwargs...)
Heat(::Val{:H}; kwargs...) = Heat(Enthalpy(); kwargs...)
Heat(::Val{:T}; kwargs...) = Heat(Temperature(); kwargs...)
Heat(para::Enthalpy; freezecurve=FreeWater(), prop=ThermalProperties(), dtlim=nothing, init=nothing) = Heat(para, prop, freezecurve, dtlim, init)
Heat(para::Temperature; freezecurve, prop=ThermalProperties(), dtlim=CFL(), init=nothing) = Heat(para, prop, freezecurve, dtlim, init)

# getter functions
thermalproperties(heat::Heat) = heat.prop
freezecurve(heat::Heat) = heat.freezecurve

# Default implementation of `variables` for freeze curve
CryoGrid.variables(::SubSurface, ::Heat, ::FreezeCurve) = ()

export HeatBC, ConstantTemp, GeothermalHeatFlux, TemperatureGradient, NFactor
include("heat_bc.jl")

export heatconduction!, enthalpy, enthalpyinv, waterice, liquidwater, freezethaw!, heatcapacity, heatcapacity!, thermalconductivity, thermalconductivity!
include("heat.jl")

end
