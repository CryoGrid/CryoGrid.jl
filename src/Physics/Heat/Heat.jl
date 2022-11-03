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
# getter functions
operator(heat::HeatBalance) = heat.op
heatproperties(heat::HeatBalance) = heat.prop
freezecurve(heat::HeatBalance) = heat.freezecurve
dtlim(heat::HeatBalance) = heat.dtlim
# default step limiters
default_dtlim(::Temperature) = Physics.CFL(maxdelta=Physics.MaxDelta(Inf))
default_dtlim(::Enthalpy) = Physics.MaxDelta(100u"kJ")
default_dtlim(::HeatOperator) = nothing

# thermal conductivity and heat capacity methods
include("thermcond.jl")
include("heatcap.jl")

# Heat operators
"""
    Diffusion{progvar,Tcond,Thc} <: HeatOperator

Represents a standard method-of-lines (MOL) forward diffusion operator for heat conduction with prognostic
`progvar`, typically either temperature `:T` or enthalpy (internal energy) `:H`.
"""
struct Diffusion{progvar,Tcond,Thc} <: HeatOperator{progvar}
    cond::Tcond
    hc::Thc
    Diffusion(progvar::Symbol, cond=quadratic_parallel_conductivity, hc=weighted_average_heatcapacity) = new{progvar,typeof(cond),typeof(hc)}(cond, hc)
end
thermalconductivity(op::Diffusion) = op.cond
heatcapacity(op::Diffusion) = op.hc
# define constructor to allow for automatic reconstruction
ConstructionBase.constructorof(::Type{<:Diffusion{progvar}}) where {progvar} = (cond,hc) -> Diffusion(progvar, cond, hc)
"""
    EnthalpyImplicit <: HeatOperator{:H}

Implicit enthalpy formulation of Swaminathan and Voller (1992) and Langer et al. (2022). Note that this
heat operator formulation does not compute a divergence `∂H∂t` but only computes the necessary diffusion
coefficients for use by an appropriate solver. See the `Solvers.LiteImplicit` module for the appropriate
solver algorithms.
"""
struct EnthalpyImplicit{Tcond,Thc} <: HeatOperator{:H}
    cond::Tcond
    hc::Thc
    EnthalpyImplicit(cond=quadratic_parallel_conductivity, hc=weighted_average_heatcapacity) = new{typeof(cond),typeof(hc)}(cond, hc)
end
thermalconductivity(op::EnthalpyImplicit) = op.cond
heatcapacity(op::EnthalpyImplicit) = op.hc

# convenience constructors for specifying prognostic variable as symbol
HeatBalance(var::Symbol=:H; kwargs...) = HeatBalance(Val{var}(); kwargs...)
HeatBalance(::Val{:H}; kwargs...) = HeatBalance(Diffusion(:H); kwargs...)
HeatBalance(::Val{:T}; kwargs...) = HeatBalance(Diffusion(:T); kwargs...)
HeatBalance(op; freezecurve=FreeWater(), prop=HeatProperties(), dtlim=default_dtlim(op)) = HeatBalance(op, prop, deepcopy(freezecurve), dtlim)
HeatBalance(op::Temperature; freezecurve, prop=HeatProperties(), dtlim=default_dtlim(op)) = HeatBalance(op, prop, deepcopy(freezecurve), dtlim)

"""
    ThermalProperties

Material thermal properties.
"""
Utils.@properties ThermalProperties(
    kh_w = 0.57u"W/m/K", # thermal conductivity of water [Hillel (1982)]
    kh_i = 2.2u"W/m/K", # thermal conductivity of ice [Hillel (1982)]
    kh_a = 0.025u"W/m/K", # thermal conductivity of air [Hillel (1982)]
    hc_w = 4.2e6u"J/K/m^3", # heat capacity of water
    hc_i = 1.9e6u"J/K/m^3", # heat capacity of ice
    hc_a = 0.00125e6u"J/K/m^3", # heat capacity of air
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

# Helper methods
"""
    enthalpy(T, C, L, θ) = T*C + L*θ

Discrete enthalpy function on temperature, heat capacity, specific latent heat of fusion, and liquid water content.
"""
@inline enthalpy(T, C, L, θ) = T*C + L*θ
"""
    enthalpyinv(H, C, L, θ) = (H - L*θ) / C

Discrete inverse enthalpy function given H, C, L, and θ.
"""
@inline enthalpyinv(H, C, L, θ) = (H - L*θ) / C
"""
    C_eff(T, C, L, ∂θw∂T, hc_w, hc_i) = C + ∂θw∂T*(L + T*(hc_w - hc_i))

Computes the apparent or "effective" heat capacity `∂H∂T` as a function of temperature, volumetric heat capacity,
latent heat of fusion, derivative of the freeze curve `∂θw∂T`, and the constituent heat capacities of water and ice.
"""
@inline C_eff(T, C, L, ∂θw∂T, hc_w, hc_i) = C + ∂θw∂T*(L + T*(hc_w - hc_i))
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

include("water_heat_coupled.jl")

end
