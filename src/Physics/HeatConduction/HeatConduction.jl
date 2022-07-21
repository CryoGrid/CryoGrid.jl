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
using FreezeCurves: FreezeCurves, FreezeCurve, FreeWater
using ModelParameters
using Unitful

import CryoGrid
import CryoGrid.Physics

export Heat, TemperatureProfile
export FreeWater, FreezeCurve, freezecurve

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

@Base.kwdef struct ThermalProperties{Tconsts,TL,Tkw,Tki,Tka,Tcw,Tci,Tca}
    consts::Tconsts = Physics.Constants()
    L::TL = consts.ρw*consts.Lsl
    kw::Tkw = 0.57u"W/m/K" # thermal conductivity of water [Hillel(1982)]
    ki::Tki = 2.2u"W/m/K" # thermal conductivity of ice [Hillel(1982)]
    ka::Tka = 0.025u"W/m/K" # air [Hillel(1982)]
    cw::Tcw = 4.2e6u"J/K/m^3" # heat capacity of water
    ci::Tci = 1.9e6u"J/K/m^3" # heat capacity of ice
    ca::Tca = 0.00125e6u"J/K/m^3" # heat capacity of air
end

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
Heat(para::Enthalpy; freezecurve=FreeWater(), prop=ThermalProperties(), dtlim=nothing, init=nothing) = Heat(para, prop, deepcopy(freezecurve), dtlim, init)
Heat(para::Temperature; freezecurve, prop=ThermalProperties(), dtlim=CFL(), init=nothing) = Heat(para, prop, deepcopy(freezecurve), dtlim, init)

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
