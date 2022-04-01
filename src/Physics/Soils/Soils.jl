module Soils

import CryoGrid: SubSurface, Parameterization
import CryoGrid: initialcondition!, variables
import CryoGrid.Physics.HeatConduction: Enthalpy, Temperature, totalwater, liquidwater, thermalconductivity, heatcapacity

using CryoGrid.Numerics
using CryoGrid.Numerics: heaviside
using CryoGrid.Physics.HeatConduction
using CryoGrid.Physics.WaterBalance
using CryoGrid.Utils

using Base: @propagate_inbounds
using IfElse
using ModelParameters
using Parameters
using Unitful

import Interpolations
import Flatten: @flattenable, flattenable

export Soil, SoilParameterization, CharacteristicFractions, SoilProfile
export soilcomponent, porosity, mineral, organic

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
abstract type SoilParameterization <: Parameterization end
"""
    CharacteristicFractions{P1,P2,P3,P4} <: SoilParameterization

Represents uniform composition of a soil volume in terms of fractions: excess ice, natural porosity, saturation, and organic solid fraction.
"""
@with_kw struct CharacteristicFractions{P1,P2,P3,P4} <: SoilParameterization
    xic::P1 = 0.0 # excess ice fraction
    por::P2 = 0.5 # natural porosity
    sat::P3 = 1.0 # saturation
    org::P4 = 0.0 # organic fraction of solid; mineral fraction is 1-org
end
# Type alias for CharacteristicFractions with all scalar/numeric constituents
const HomogeneousCharacteristicFractions = CharacteristicFractions{<:Number,<:Number,<:Number,<:Number}
SoilProfile(pairs::Pair{<:DistQuantity,<:SoilParameterization}...) = Profile(pairs...)
# Helper functions for obtaining soil compositions from characteristic fractions.
soilcomponent(::Val{var}, para::CharacteristicFractions) where var = soilcomponent(Val{var}(), para.xic, para.por, para.sat, para.org)
soilcomponent(::Val{:θp}, χ, ϕ, θ, ω) = (1-χ)*ϕ
soilcomponent(::Val{:θw}, χ, ϕ, θ, ω) = χ + (1-χ)*ϕ*θ
soilcomponent(::Val{:θm}, χ, ϕ, θ, ω) = (1-χ)*(1-ϕ)*(1-ω)
soilcomponent(::Val{:θo}, χ, ϕ, θ, ω) = (1-χ)*(1-ϕ)*ω
"""
Thermal conductivity constants.
"""
@with_kw struct SoilTCParams @deftype Float"W/(m*K)"
    kw = 0.57xu"W/(m*K)" #water [Hillel(1982)]
    ko = 0.25xu"W/(m*K)" #organic [Hillel(1982)]
    km = 3.8xu"W/(m*K)" #mineral [Hillel(1982)]
    ka = 0.025xu"W/(m*K)" #air [Hillel(1982)]
    ki = 2.2xu"W/(m*K)" #ice [Hillel(1982)]
end
"""
Heat capacity constants.
"""
@with_kw struct SoilHCParams @deftype Float"J/(K*m^3)"
    cw = 4.2e6xu"J/(K*m^3)" #[J/m^3K] heat capacity water
    co = 2.5e6xu"J/(K*m^3)" #[J/m^3K]  heat capacity organic
    cm = 2e6xu"J/(K*m^3)" #[J/m^3K]  heat capacity mineral
    ca = 0.00125e6xu"J/(K*m^3)" #[J/m^3K]  heat capacity pore space
    ci = 1.9e6xu"J/(K*m^3)" #[J/m^3K]  heat capacity ice
end
"""
Basic Soil layer.
"""
@with_kw struct Soil{P<:SoilParameterization,S} <: SubSurface
    para::P = CharacteristicFractions()
    tc::SoilTCParams = SoilTCParams()
    hc::SoilHCParams = SoilHCParams()
    sp::S = nothing # user-defined specialization
end
# SoilComposition trait impl
SoilComposition(soil::Soil) = SoilComposition(typeof(soil))
SoilComposition(::Type{<:Soil}) = Heterogeneous()
SoilComposition(::Type{<:Soil{<:HomogeneousCharacteristicFractions}}) = Homogeneous()
# Fixed volumetric content methods
totalwater(soil::Soil, state) = totalwater(SoilComposition(soil), soil, state)
porosity(soil::Soil, state) = porosity(SoilComposition(soil), soil, state)
mineral(soil::Soil, state) = mineral(SoilComposition(soil), soil, state)
organic(soil::Soil, state) = organic(SoilComposition(soil), soil, state)
## Homogeneous soils
totalwater(::Homogeneous, soil::Soil{<:CharacteristicFractions}, state) = soilcomponent(Val{:θw}(), soil.para)
porosity(::Homogeneous, soil::Soil{<:CharacteristicFractions}, state) = soilcomponent(Val{:θp}(), soil.para)
mineral(::Homogeneous, soil::Soil{<:CharacteristicFractions}, state) = soilcomponent(Val{:θm}(), soil.para)
organic(::Homogeneous, soil::Soil{<:CharacteristicFractions}, state) = soilcomponent(Val{:θo}(), soil.para)
## Heterogeneous soils
totalwater(::Heterogeneous, soil::Soil, state) = state.θw   
porosity(::Heterogeneous, soil::Soil, state) = state.θp
mineral(::Heterogeneous, soil::Soil, state) = state.θm
organic(::Heterogeneous, soil::Soil, state) = state.θo

variables(soil::Soil) = variables(SoilComposition(soil), soil)
variables(::Homogeneous, ::Soil) = ()
variables(::Heterogeneous, ::Soil) = (
    Diagnostic(:θw, Float64, OnGrid(Cells)),
    Diagnostic(:θp, Float64, OnGrid(Cells)),
    Diagnostic(:θm, Float64, OnGrid(Cells)),
    Diagnostic(:θo, Float64, OnGrid(Cells)),
)

initialcondition!(soil::Soil, state) = initialcondition!(SoilComposition(soil), soil, state)
initialcondition!(::Homogeneous, ::Soil, state) = nothing
"""
    initialcondition!(::Heterogeneous, soil::Soil{<:CharacteristicFractions}, state)

Default implementation of initialcondition! for heterogeneous soils parameterized by characteristic fractions.
Fields of `soil.para` may be either numbers (including `Param`s) or `Function`s of the form `f(::Soil, state)::AbstractVector`.
"""
function initialcondition!(::Heterogeneous, soil::Soil{<:CharacteristicFractions}, state)
    evaluate(x::Number) = x
    evaluate(f::Function) = f(soil, state)
    χ = evaluate(soil.para.xic)
    ϕ = evaluate(soil.para.por)
    θ = evaluate(soil.para.sat)
    ω = evaluate(soil.para.org)
    @. state.θw = soilcomponent(Val{:θw}(), χ, ϕ, θ, ω)
    @. state.θp = soilcomponent(Val{:θp}(), χ, ϕ, θ, ω)
    @. state.θm = soilcomponent(Val{:θm}(), χ, ϕ, θ, ω)
    @. state.θo = soilcomponent(Val{:θo}(), χ, ϕ, θ, ω)
end

export SFCC, DallAmico, Westermann, McKenzie, SFCCNewtonSolver, SFCCPreSolver
include("soilheat.jl")

end
