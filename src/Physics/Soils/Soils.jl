module Soils

using CryoGrid: SubSurface, Parameterization
using CryoGrid.Numerics
using CryoGrid.Numerics: heaviside
using CryoGrid.Physics
using CryoGrid.Physics.HeatConduction
using CryoGrid.Physics.Hydrology
using CryoGrid.Utils

using Base: @propagate_inbounds, @kwdef
using IfElse
using Interpolations: Interpolations
using IntervalSets
using ForwardDiff
using ModelParameters
using Unitful

import CryoGrid
import CryoGrid.Physics
import CryoGrid.Physics.HeatConduction

import Flatten: @flattenable, flattenable

export Soil, SoilParameterization, CharacteristicFractions, SoilProfile
export soilparameters, soilcomponent, porosity, mineral, organic

const Enthalpy = HeatConduction.Enthalpy
const Temperature = HeatConduction.Temperature

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
    xic::P1 = Param(0.0, domain=0..1) # excess ice fraction
    por::P2 = Param(0.5, domain=0..1) # natural porosity
    sat::P3 = Param(1.0, domain=0..1) # saturation
    org::P4 = Param(0.0, domain=0..1) # organic fraction of solid; mineral fraction is 1-org
end
soilparameters(::Type{CharacteristicFractions}=CharacteristicFractions; xic, por, sat, org) = CharacteristicFractions(xic=Param(xic, domain=0..1), por=Param(por, domain=0..1), sat=Param(sat, domain=0..1), org=Param(org, domain=0..1))
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
SoilThermalProperties(;
    ko = Param(0.25, units=u"W/m/K"), # organic [Hillel(1982)]
    km = Param(3.8, units=u"W/m/K"), # mineral [Hillel(1982)]
    co = Param(2.5e6, units=u"J/K/m^3"), # heat capacity organic
    cm = Param(2.0e6, units=u"J/K/m^3"), # heat capacity mineral
) = (; ko, km, co, cm)
"""
Basic Soil layer.
"""
@kwdef struct Soil{TPara<:SoilParameterization,TProp,TSp} <: SubSurface
    para::TPara = CharacteristicFractions()
    prop::TProp = SoilThermalProperties()
    sp::TSp = nothing # user-defined specialization
end
HeatConduction.thermalproperties(soil::Soil) = soil.prop
# SoilComposition trait impl
SoilComposition(soil::Soil) = SoilComposition(typeof(soil))
SoilComposition(::Type{<:Soil}) = Heterogeneous()
SoilComposition(::Type{<:Soil{<:HomogeneousCharacteristicFractions}}) = Homogeneous()
# Volumetric fraction methods
Physics.waterice(soil::Soil, state) = waterice(SoilComposition(soil), soil, state)
porosity(soil::Soil, state) = porosity(SoilComposition(soil), soil, state)
mineral(soil::Soil, state) = mineral(SoilComposition(soil), soil, state)
organic(soil::Soil, state) = organic(SoilComposition(soil), soil, state)
## Homogeneous soils
Physics.waterice(::Homogeneous, soil::Soil{<:CharacteristicFractions}, state=nothing) = soilcomponent(Val{:θwi}(), soil.para)
porosity(::Homogeneous, soil::Soil{<:CharacteristicFractions}, state=nothing) = soilcomponent(Val{:θp}(), soil.para)
mineral(::Homogeneous, soil::Soil{<:CharacteristicFractions}, state=nothing) = soilcomponent(Val{:θm}(), soil.para)
organic(::Homogeneous, soil::Soil{<:CharacteristicFractions}, state=nothing) = soilcomponent(Val{:θo}(), soil.para)
## Heterogeneous soils
Physics.waterice(::Heterogeneous, soil::Soil, state) = state.θwi   
porosity(::Heterogeneous, soil::Soil, state) = state.θp
mineral(::Heterogeneous, soil::Soil, state) = state.θm
organic(::Heterogeneous, soil::Soil, state) = state.θo

CryoGrid.variables(soil::Soil) = CryoGrid.variables(SoilComposition(soil), soil)
CryoGrid.variables(::Homogeneous, ::Soil) = ()
CryoGrid.variables(::Heterogeneous, ::Soil) = (
    Diagnostic(:θwi, OnGrid(Cells), domain=0..1),
    Diagnostic(:θp, OnGrid(Cells), domain=0..1),
    Diagnostic(:θm, OnGrid(Cells), domain=0..1),
    Diagnostic(:θo, OnGrid(Cells), domain=0..1),
)

CryoGrid.initialcondition!(soil::Soil, state) = CryoGrid.initialcondition!(SoilComposition(soil), soil, state)
CryoGrid.initialcondition!(::Homogeneous, ::Soil, state) = nothing
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

export SWRC, VanGenuchten
include("sfcc/swrc.jl")
export SFCC, DallAmico, Westermann, McKenzie, SFCCNewtonSolver, SFCCPreSolver
include("sfcc/sfcc.jl")
include("soilheat.jl")

end
