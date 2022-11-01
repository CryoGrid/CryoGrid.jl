module Heat

import CryoGrid
import CryoGrid.Physics

import ConstructionBase

using CryoGrid
using CryoGrid.InputOutput: Forcing
using CryoGrid.Physics
using CryoGrid.Physics.Boundaries
using CryoGrid.Physics.Hydrology
using CryoGrid.Numerics
using CryoGrid.Numerics: nonlineardiffusion!, harmonicmean!, harmonicmean, heaviside
using CryoGrid.Utils

using Base: @propagate_inbounds
using IfElse
using FreezeCurves: FreezeCurves, FreezeCurve, FreeWater
using ModelParameters
using Unitful
using UnPack

export HeatBalance, HeatProperties, ThermalProperties, TemperatureProfile
export FreeWater, FreezeCurve
export thermalproperties, freezecurve

"""
    HeatOperator{progvar}

Base type for different numerical formulations of heat conduction.
"""
abstract type HeatOperator{progvar} end
"""
    Diffusion{progvar} <: HeatOperator

Represents a standard method-of-lines (MOL) forward diffusion operator for heat conduction with prognostic
`progvar`, typically either temperature `:T` or enthalpy (internal energy) `:H`.
"""
struct Diffusion{progvar} <: HeatOperator{progvar}
    Diffusion(progvar::Symbol) = new{progvar}()
end
# define constructor to allow for automatic reconstruction
ConstructionBase.constructorof(::Type{Diffusion{progvar}}) where {progvar} = () -> Diffusion(progvar)
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
const Enthalpy = HeatOperator{:H}
"""
Type alias for `HeatOperator{:T}`, i.e. temperature-based heat conduction operators.
"""
const Temperature = HeatOperator{:T}
"""
Numerical constants for pararameterizing heat processes.
"""
Utils.@properties HeatProperties(
    ρw = Physics.Constants.ρw,
    Lsl = Physics.Constants.Lsl,
    L = ρw*Lsl,
)
# do not parameterize heat properties
CryoGrid.parameterize(prop::HeatProperties) = prop
"""
    HeatBalance{Tfc<:FreezeCurve,THeatOp<:HeatOperator,Tdt,Tprop} <: SubSurfaceProcess

Represents subsurface heat transfer processes. The formulation of heat transfer is governed by
the `HeatOperator`, `op` and 
"""
struct HeatBalance{Tfc<:FreezeCurve,THeatOp<:HeatOperator,Tdt,Tprop} <: SubSurfaceProcess
    op::THeatOp
    prop::Tprop
    freezecurve::Tfc
    dtlim::Tdt  # timestep limiter
end
default_dtlim(::Temperature) = Physics.CFL(maxdelta=Physics.MaxDelta(Inf))
default_dtlim(::Enthalpy) = Physics.MaxDelta(100u"kJ")
default_dtlim(::HeatOperator) = nothing
# convenience constructors for specifying prognostic variable as symbol
HeatBalance(var::Symbol=:H; kwargs...) = HeatBalance(Val{var}(); kwargs...)
HeatBalance(::Val{:H}; kwargs...) = HeatBalance(Diffusion(:H); kwargs...)
HeatBalance(::Val{:T}; kwargs...) = HeatBalance(Diffusion(:T); kwargs...)
HeatBalance(op; freezecurve=FreeWater(), prop=HeatProperties(), dtlim=default_dtlim(op)) = HeatBalance(op, prop, deepcopy(freezecurve), dtlim)
HeatBalance(op::Temperature; freezecurve, prop=HeatProperties(), dtlim=default_dtlim(op)) = HeatBalance(op, prop, deepcopy(freezecurve), dtlim)

# getter functions
freezecurve(heat::HeatBalance) = heat.freezecurve

"""
    ThermalProperties

Material thermal properties.
"""
Utils.@properties ThermalProperties(
    kw = 0.57u"W/m/K", # thermal conductivity of water [Hillel (1982)]
    ki = 2.2u"W/m/K", # thermal conductivity of ice [Hillel (1982)]
    ka = 0.025u"W/m/K", # thermal conductivity of air [Hillel (1982)]
    cw = 4.2e6u"J/K/m^3", # heat capacity of water
    ci = 1.9e6u"J/K/m^3", # heat capacity of ice
    ca = 0.00125e6u"J/K/m^3", # heat capacity of air
)
function CryoGrid.parameterize(prop::ThermalProperties)
    return ThermalProperties(
        map(values(prop)) do val
            # this currently assumes that all properties have a strictly positive domain!
            CryoGrid.parameterize(val, domain=StrictlyPositive)
        end
    )
end
thermalproperties(::TLayer) where {TLayer<:SubSurface} = error("thermal properties not defined for $TLayer")

"""
    TemperatureProfile(pairs::Pair{<:Union{DistQuantity,Param},<:Union{TempQuantity,Param}}...)

Convenience constructor for `Numerics.Profile` which automatically converts temperature quantities.
"""
TemperatureProfile(pairs::Pair{<:Union{DistQuantity,Param},<:Union{TempQuantity,Param}}...) = Profile(map(p -> uconvert(u"m", p[1]) => uconvert(u"°C", p[2]), pairs))

export HeatBC, ConstantTemperature, GeothermalHeatFlux, TemperatureGradient, NFactor
include("heat_bc.jl")

export heatconduction!, enthalpy, enthalpyinv, freezethaw!, heatcapacity, heatcapacity!, thermalconductivity, thermalconductivity!
include("heat_conduction.jl")

export ImplicitHeat
include("heat_implicit.jl")

include("water_heat.jl")

end
