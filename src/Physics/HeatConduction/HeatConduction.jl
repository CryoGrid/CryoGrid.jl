module HeatConduction

import CryoGrid
import CryoGrid.Physics

using CryoGrid
using CryoGrid.InputOutput: Forcing
using CryoGrid.Physics
using CryoGrid.Physics.Boundaries
using CryoGrid.Physics.Hydrology
using CryoGrid.Numerics
using CryoGrid.Numerics: nonlineardiffusion!, harmonicmean!, harmonicmean, heaviside
using CryoGrid.Utils

using Base: @propagate_inbounds, @kwdef
using IfElse
using FreezeCurves: FreezeCurves, FreezeCurve, FreeWater
using ModelParameters
using Unitful

export Heat, TemperatureProfile
export FreeWater, FreezeCurve, freezecurve

"""
    TemperatureProfile(pairs::Pair{<:Union{DistQuantity,Param},<:Union{TempQuantity,Param}}...)

Convenience constructor for `Numerics.Profile` which automatically converts temperature quantities.
"""
TemperatureProfile(pairs::Pair{<:Union{DistQuantity,Param},<:Union{TempQuantity,Param}}...) = Profile(map(p -> uconvert(u"m", p[1]) => uconvert(u"°C", p[2]), pairs))

"""
    HeatFormulation

Base type for different numerical formulations of two-phase heat diffusion.
"""
abstract type HeatFormulation end
struct Enthalpy <: HeatFormulation end
struct Temperature <: HeatFormulation end

@Base.kwdef struct ThermalProperties{Tconsts,TL,Tkw,Tki,Tka,Tcw,Tci,Tca}
    consts::Tconsts = Physics.Constants()
    L::TL = consts.ρw*consts.Lsl
    kw::Tkw = Param(0.57, units=u"W/m/K", domain=StrictlyPositive) # thermal conductivity of water [Hillel (1982)]
    ki::Tki = Param(2.2, units=u"W/m/K", domain=StrictlyPositive) # thermal conductivity of ice [Hillel (1982)]
    ka::Tka = Param(0.025, unit=u"W/m/K", domain=StrictlyPositive) # thermal conductivity of air [Hillel (1982)]
    cw::Tcw = Param(4.2e6, units=u"J/K/m^3", domain=StrictlyPositive) # heat capacity of water
    ci::Tci = Param(1.9e6, units=u"J/K/m^3", domain=StrictlyPositive) # heat capacity of ice
    ca::Tca = Param(0.00125e6, units=u"J/K/m^3", domain=StrictlyPositive) # heat capacity of air
end

struct Heat{Tfc<:FreezeCurve,TForm<:HeatFormulation,Tdt,Tinit,TProp} <: SubSurfaceProcess
    form::TForm
    prop::TProp
    freezecurve::Tfc
    dtlim::Tdt  # timestep limiter
    init::Tinit # optional initialization scheme
end
const DEFAULT_MAX_ENERGY_CHANGE = 50u"kJ"
_default_dtlim(::Union{Enthalpy,Temperature}, ::FreezeCurve) = Physics.CFL(maxdelta=Physics.MaxDelta(DEFAULT_MAX_ENERGY_CHANGE))
_default_dtlim(::Enthalpy, ::FreeWater) = Physics.MaxDelta(DEFAULT_MAX_ENERGY_CHANGE) # CFL doesn't work with FreeWater freeze curve
# convenience constructors for specifying prognostic variable as symbol
Heat(var::Symbol=:H; kwargs...) = Heat(Val{var}(); kwargs...)
Heat(::Val{:H}; kwargs...) = Heat(Enthalpy(); kwargs...)
Heat(::Val{:T}; kwargs...) = Heat(Temperature(); kwargs...)
Heat(form::Enthalpy; freezecurve=FreeWater(), prop=ThermalProperties(), dtlim=_default_dtlim(form, freezecurve), init=nothing) = Heat(form, prop, deepcopy(freezecurve), dtlim, init)
Heat(form::Temperature; freezecurve, prop=ThermalProperties(), dtlim=_default_dtlim(form, freezecurve), init=nothing) = Heat(form, prop, deepcopy(freezecurve), dtlim, init)

# getter functions
thermalproperties(heat::Heat) = heat.prop
freezecurve(heat::Heat) = heat.freezecurve

export HeatBC, ConstantTemp, GeothermalHeatFlux, TemperatureGradient, NFactor
include("heat_bc.jl")

export heatconduction!, enthalpy, enthalpyinv, freezethaw!, heatcapacity, heatcapacity!, thermalconductivity, thermalconductivity!
include("heat.jl")

end
