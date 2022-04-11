module HeatConduction

import CryoGrid: SubSurfaceProcess, BoundaryStyle, Dirichlet, Neumann, BoundaryProcess, Layer, Top, Bottom, SubSurface, Callback
import CryoGrid: diagnosticstep!, prognosticstep!, interact!, initialcondition!, boundaryflux, boundaryvalue, variables, callbacks, criterion, affect!

using CryoGrid.InputOutput: Forcing
using CryoGrid.Physics
using CryoGrid.Physics.Boundaries
using CryoGrid.Numerics
using CryoGrid.Numerics: nonlineardiffusion!, harmonicmean!, harmonicmean, heaviside
using CryoGrid.Utils

using Base: @propagate_inbounds, @kwdef
using IfElse
using ModelParameters
using SimulationLogs
using Unitful

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
    HeatImpl

Base type for different numerical formulations of two-phase heat diffusion.
"""
abstract type HeatImpl end
struct Enthalpy <: HeatImpl end
struct Temperature <: HeatImpl end

"""
    ThermalProperties

Base type for defining thermal properties.
"""
abstract type ThermalProperties <: IterableStruct end
"""
    HydroThermalProperties{Tρ,TLf,Tkw,Tki,Tcw,Tci}

Thermal properties of water used in two-phase heat conduction.
"""
@kwdef struct HydroThermalProperties{Tρ,TLf,Tkw,Tki,Tcw,Tci} <: ThermalProperties
    ρ::Tρ = Param(1000.0, units=u"kg/m^3") # density of water
    Lf::TLf = Param(334000.0, units=u"J/kg") # latent heat of fusion of water
    kw::Tkw = Param(0.57, units=u"W/m/K") # thermal conductivity of water [Hillel(1982)]
    ki::Tki = Param(2.2, units=u"W/m/K") # thermal conductivity of ice [Hillel(1982)]
    cw::Tcw = Param(4.2e6, units=u"J/K/m^3") # heat capacity of water
    ci::Tci = Param(1.9e6, units=u"J/K/m^3") # heat capacity of ice
end

@kwdef struct Heat{Tfc<:FreezeCurve,Tsp,Tinit,TProp<:HydroThermalProperties,TL} <: SubSurfaceProcess
    prop::TProp = HydroThermalProperties()
    L::TL = prop.ρ*prop.Lf # [J/m^3] (specific latent heat of fusion of water)
    freezecurve::Tfc = FreeWater() # freeze curve, defautls to free water fc
    sp::Tsp = Enthalpy() # specialization
    init::Tinit = nothing # optional initialization scheme
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
