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
    HeatOperator{progvar}

Base type for different numerical formulations of heat conduction.
"""
abstract type HeatOperator{progvar} end
"""
    Diffusion{progvar} <: HeatOperator

Standard method-of-lines (MOL) forward diffusion operator for heat conduction with prognostic
`progvar`, typically either temperature `:T` or enthalpy (internal energy) `:H`.
"""
struct Diffusion{progvar} <: HeatOperator{progvar}
    Diffusion(progvar::Symbol) = new{progvar}()
end
"""
    EnthalpyImplicit <: HeatOperator{:H}

Implicit enthalpy formulation of Swaminathan and Voller (1992) and Langer et al. (2022). Note that this
heat operator formulation does not compute a divergence `∂H∂t` but only computes the necessary diffusion
coefficients for use by an appropriate solver. See the `Solvers.LiteImplicit` module for the appropriate
solver algorithms.
"""
struct EnthalpyImplicit <: HeatOperator{:H} end
"""
Type alias for `HeatOperator{:H}`, i.e. enthalpy-based heat conduction operators.
"""
const PrognosticEnthalpy = HeatOperator{:H}
"""
Type alias for `HeatOperator{:T}`, i.e. temperature-based heat conduction operators.
"""
const PrognosticTemperature = HeatOperator{:T}

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

struct Heat{Tfc<:FreezeCurve,THeatOp<:HeatOperator,Tdt,Tinit,TProp} <: SubSurfaceProcess
    op::THeatOp
    prop::TProp
    freezecurve::Tfc
    dtlim::Tdt  # timestep limiter
    init::Tinit # optional initialization scheme
end
default_dtlim(::PrognosticTemperature, ::FreezeCurve) = Physics.CFL(maxdelta=Physics.MaxDelta(Inf))
default_dtlim(::PrognosticEnthalpy, ::FreezeCurve) = Physics.MaxDelta(100u"kJ")
default_dtlim(::HeatOperator, ::FreezeCurve) = nothing
# convenience constructors for specifying prognostic variable as symbol
Heat(var::Symbol=:H; kwargs...) = Heat(Val{var}(); kwargs...)
Heat(::Val{:H}; kwargs...) = Heat(Diffusion(:H); kwargs...)
Heat(::Val{:T}; kwargs...) = Heat(Diffusion(:T); kwargs...)
Heat(op; freezecurve=FreeWater(), prop=ThermalProperties(), dtlim=default_dtlim(op, freezecurve), init=nothing) = Heat(op, prop, deepcopy(freezecurve), dtlim, init)
Heat(op::PrognosticTemperature; freezecurve, prop=ThermalProperties(), dtlim=default_dtlim(op, freezecurve), init=nothing) = Heat(op, prop, deepcopy(freezecurve), dtlim, init)

# getter functions
thermalproperties(heat::Heat) = heat.prop
freezecurve(heat::Heat) = heat.freezecurve

"""
    TemperatureProfile(pairs::Pair{<:Union{DistQuantity,Param},<:Union{TempQuantity,Param}}...)

Convenience constructor for `Numerics.Profile` which automatically converts temperature quantities.
"""
TemperatureProfile(pairs::Pair{<:Union{DistQuantity,Param},<:Union{TempQuantity,Param}}...) = Profile(map(p -> uconvert(u"m", p[1]) => uconvert(u"°C", p[2]), pairs))

export HeatBC, ConstantTemperature, GeothermalHeatFlux, TemperatureGradient, NFactor
include("heat_bc.jl")

export heatconduction!, enthalpy, enthalpyinv, freezethaw!, heatcapacity, heatcapacity!, thermalconductivity, thermalconductivity!
include("heat.jl")

export ImplicitHeat
include("heat_implicit.jl")

include("water_heat.jl")

end
