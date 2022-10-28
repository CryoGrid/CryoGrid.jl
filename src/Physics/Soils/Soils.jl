module Soils

using CryoGrid
using CryoGrid.Numerics
using CryoGrid.Physics
using CryoGrid.Physics.HeatConduction
using CryoGrid.Physics.Hydrology
using CryoGrid.Utils

using Base: @propagate_inbounds, @kwdef
using IfElse
using IntervalSets
using ForwardDiff
using FreezeCurves
using ModelParameters
using Setfield
using StaticArrays
using Unitful

import CryoGrid
import CryoGrid.InputOutput
import CryoGrid.Physics
import CryoGrid.Physics.HeatConduction
import CryoGrid.Physics.Hydrology

export Soil, SoilParameterization, CharacteristicFractions, SoilProfile
export soilparameters, soilcomponent, porosity, mineral, organic

# from FreezeCurves
export SFCC, PainterKarra, DallAmico, DallAmicoSalt, Westermann, McKenzie, VanGenuchten, BrooksCorey

# aliases for heat formulations in HeatConduction module
const Temperature = HeatConduction.Temperature
const Enthalpy = HeatConduction.Enthalpy
const EnthalpyImplicit = HeatConduction.EnthalpyImplicit

"""
    SoilComposition

Trait for representing homogenous vs heterogeneous soil layers.
"""
abstract type SoilComposition end
struct Homogeneous <: SoilComposition end
struct Heterogeneous <: SoilComposition end
"""
    SoilParameterization

Abstract base type for parameterizations of soil properties.
"""
abstract type SoilParameterization end
"""
    CharacteristicFractions{P1,P2,P3,P4} <: SoilParameterization

Represents uniform composition of a soil volume in terms of fractions: excess ice, natural porosity, saturation, and organic solid fraction.
"""
Base.@kwdef struct CharacteristicFractions{P1,P2,P3,P4} <: SoilParameterization
    xic::P1 = 0.0 # excess ice fraction
    por::P2 = 0.5 # natural porosity
    sat::P3 = 1.0 # saturation
    org::P4 = 0.0 # organic fraction of solid; mineral fraction is 1-org
end
function soilparameters(::Type{CharacteristicFractions}=CharacteristicFractions; xic, por, sat, org)
    CharacteristicFractions(
        xic=Param(xic, domain=0..1),
        por=Param(por, domain=0..1),
        sat=Param(sat, domain=0..1),
        org=Param(org, domain=0..1)
    )
end
# Type alias for CharacteristicFractions with all scalar/numeric constituents
const HomogeneousCharacteristicFractions = CharacteristicFractions{<:Number,<:Number,<:Number,<:Number}
SoilProfile(pairs::Pair{<:DistQuantity,<:SoilParameterization}...) = Profile(pairs...)
# Helper functions for obtaining soil compositions from characteristic fractions.
soilcomponent(::Val{var}, para::CharacteristicFractions) where var = soilcomponent(Val{var}(), para.xic, para.por, para.sat, para.org)
soilcomponent(::Val{:θp}, χ, ϕ, θ, ω) = (1-χ)*ϕ
soilcomponent(::Val{:θwi}, χ, ϕ, θ, ω) = χ + (1-χ)*ϕ*θ
soilcomponent(::Val{:θm}, χ, ϕ, θ, ω) = (1-χ)*(1-ϕ)*(1-ω)
soilcomponent(::Val{:θo}, χ, ϕ, θ, ω) = (1-χ)*(1-ϕ)*ω
"""
Soil thermal properties.
"""
Base.@kwdef struct SoilThermalProperties{Tko,Tkm,Tco,Tcm}
    ko::Tko = Param(0.25, units=u"W/m/K", domain=StrictlyPositive) # organic [Hillel (1982)]
    km::Tkm = Param(3.8, units=u"W/m/K", domain=StrictlyPositive) # mineral [Hillel (1982)]
    co::Tco = Param(2.5e6, units=u"J/K/m^3", domain=StrictlyPositive) # heat capacity organic
    cm::Tcm = Param(2.0e6, units=u"J/K/m^3", domain=StrictlyPositive) # heat capacity mineral
end
"""
    Soil{Tpara<:SoilParameterization,Tprop,Tsp,Tproc} <: SubSurface{Tproc}

Generic Soil layer.
"""
@kwdef struct Soil{Tpara<:SoilParameterization,Tprop,Tsp,Tproc} <: SubSurface{Tproc}
    para::Tpara = CharacteristicFractions()
    prop::Tprop = SoilThermalProperties()
    sp::Tsp = nothing # user-defined specialization
    proc::Tproc
end
Soil(proc::Tproc; kwargs...) where Tproc = Soil(;proc, kwargs...)
HeatConduction.thermalproperties(soil::Soil) = soil.prop
# SoilComposition trait impl
SoilComposition(soil::Soil) = SoilComposition(typeof(soil))
SoilComposition(::Type{<:Soil}) = Heterogeneous()
SoilComposition(::Type{<:Soil{<:HomogeneousCharacteristicFractions}}) = Homogeneous()
# Volumetric fraction methods
porosity(soil::Soil, state) = porosity(SoilComposition(soil), soil, state)
mineral(soil::Soil, state) = mineral(SoilComposition(soil), soil, state)
organic(soil::Soil, state) = organic(SoilComposition(soil), soil, state)
## Homogeneous soils
porosity(::Homogeneous, soil::Soil{<:CharacteristicFractions}, state=nothing) = soilcomponent(Val{:θp}(), soil.para)
mineral(::Homogeneous, soil::Soil{<:CharacteristicFractions}, state=nothing) = soilcomponent(Val{:θm}(), soil.para)
organic(::Homogeneous, soil::Soil{<:CharacteristicFractions}, state=nothing) = soilcomponent(Val{:θo}(), soil.para)
## Heterogeneous soils  
porosity(::Heterogeneous, soil::Soil, state) = state.θp
mineral(::Heterogeneous, soil::Soil, state) = state.θm
organic(::Heterogeneous, soil::Soil, state) = state.θo

CryoGrid.variables(soil::Soil) = CryoGrid.variables(SoilComposition(soil), soil)
CryoGrid.variables(::Homogeneous, soil::Soil) = (
    Diagnostic(:θwi, OnGrid(Cells), domain=0..1),
    CryoGrid.variables(soil, processes(soil))...
)
CryoGrid.variables(::Heterogeneous, soil::Soil) = (
    Diagnostic(:θwi, OnGrid(Cells), domain=0..1),
    Diagnostic(:θp, OnGrid(Cells), domain=0..1),
    Diagnostic(:θm, OnGrid(Cells), domain=0..1),
    Diagnostic(:θo, OnGrid(Cells), domain=0..1),
    CryoGrid.variables(soil, processes(soil))...,
)

CryoGrid.initialcondition!(soil::Soil, state) = CryoGrid.initialcondition!(SoilComposition(soil), soil, state)
function CryoGrid.initialcondition!(::Homogeneous, soil::Soil{<:CharacteristicFractions}, state)
    χ = soil.para.xic
    ϕ = soil.para.por
    θ = soil.para.sat
    ω = soil.para.org
    @. state.θwi = soilcomponent(Val{:θwi}(), χ, ϕ, θ, ω)
    CryoGrid.initialcondition!(soil, processes(soil), state)
end
"""
    initialcondition!(::Heterogeneous, soil::Soil{<:CharacteristicFractions}, state)

Default implementation of initialcondition! for heterogeneous soils parameterized by characteristic fractions.
Fields of `soil.para` may be either numbers (including `Param`s) or `Function`s of the form `f(::Soil, state)::AbstractVector`.
"""
function CryoGrid.initialcondition!(::Heterogeneous, soil::Soil{<:CharacteristicFractions}, state)
    evaluate(x::Number) = x
    evaluate(f::Function) = f(soil, state)
    χ = evaluate(soil.para.xic)
    ϕ = evaluate(soil.para.por)
    θ = evaluate(soil.para.sat)
    ω = evaluate(soil.para.org)
    @. state.θwi = soilcomponent(Val{:θwi}(), χ, ϕ, θ, ω)
    @. state.θp = soilcomponent(Val{:θp}(), χ, ϕ, θ, ω)
    @. state.θm = soilcomponent(Val{:θm}(), χ, ϕ, θ, ω)
    @. state.θo = soilcomponent(Val{:θo}(), χ, ϕ, θ, ω)
    CryoGrid.initialcondition!(soil, processes(soil), state)
end

"""
    mineral(soil::Soil, state, i)

Retrieves the mineral content for the given layer at grid cell `i`, if provided.
Defaults to using the scalar mineral content defined on `soil`.
"""
@inline mineral(soil::Soil, state, i) = Utils.getscalar(mineral(soil, state), i)
"""
    organic(soil::Soil, state, i)

Retrieves the organic content for the given layer at grid cell `i`, if provided.
Defaults to using the scalar organic content defined on `soil`.
"""
@inline organic(soil::Soil, state, i) = Utils.getscalar(organic(soil, state), i)
"""
    porosity(soil::Soil, state, i)

Retrieves the porosity for the given layer at grid cell `i`, if provided.
Defaults to using the scalar porosity defined on `soil`.
"""
@inline porosity(soil::Soil, state, i) = Utils.getscalar(porosity(soil, state), i)

export RichardsEq
include("soilwater.jl")
include("soilheat.jl")

end
