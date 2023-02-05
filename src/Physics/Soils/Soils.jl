module Soils

using CryoGrid
using CryoGrid.Numerics
using CryoGrid.Physics
using CryoGrid.Physics.Heat
using CryoGrid.Physics.Hydrology
using CryoGrid.Utils

using IfElse
using IntervalSets
using ForwardDiff
using FreezeCurves
using ModelParameters
using Setfield
using StaticArrays
using Unitful
using UnPack

import CryoGrid
import CryoGrid.InputOutput
import CryoGrid.Physics
import CryoGrid.Physics.Heat
import CryoGrid.Physics.Hydrology

# from FreezeCurves
export SFCC, PainterKarra, DallAmico, DallAmicoSalt, Westermann, McKenzie
export VanGenuchten, BrooksCorey

# aliases for heat formulations in Heat module
const Temperature = Heat.Temperature
const Enthalpy = Heat.Enthalpy
const EnthalpyImplicit = Heat.EnthalpyImplicit

export Soil, SoilParameterization, SoilProperties
include("types.jl")

export porosity, mineral, organic
include("methods.jl")

export HomogeneousMixture, SoilProfile, soilcomponent
include("soil_para.jl")

## Homogeneous soils
porosity(soil::Soil{<:HomogeneousMixture}, state=nothing) = soilcomponent(Val{:θp}(), soil.para)
mineral(soil::Soil{<:HomogeneousMixture}, state=nothing) = soilcomponent(Val{:θm}(), soil.para)
organic(soil::Soil{<:HomogeneousMixture}, state=nothing) = soilcomponent(Val{:θo}(), soil.para)

"""
    SoilProfile(pairs::Pair{<:DistQuantity,<:SoilParameterization}...)

Alias for `Profile(pairs...)` assigning soil parameterizations to specific depths.
"""
SoilProfile(pairs::Pair{<:DistQuantity,<:SoilParameterization}...) = Profile(pairs...)

export RichardsEq
include("soil_water.jl")
include("soil_heat.jl")
include("soil_water_heat_coupled.jl")

CryoGrid.variables(soil::Soil{<:HomogeneousMixture}) = CryoGrid.variables(soil, processes(soil))
CryoGrid.variables(soil::Soil) = (
    Diagnostic(:θwi, OnGrid(Cells), domain=0..1),
    Diagnostic(:θp, OnGrid(Cells), domain=0..1),
    Diagnostic(:θm, OnGrid(Cells), domain=0..1),
    Diagnostic(:θo, OnGrid(Cells), domain=0..1),
    CryoGrid.variables(soil, processes(soil))...,
)
function CryoGrid.initialcondition!(soil::Soil{<:MineralSediment}, state)
    evaluate(x::Number) = x
    evaluate(f::Function) = f(soil, state)
    ϕ = evaluate(soil.para.por)
    θ = evaluate(soil.para.sat)
    @. state.θp = ϕ
    @. state.θwi = ϕ*θ
    @. state.θm = 1 - state.θp
    @. state.θo = zero(eltype(state.θo))
    CryoGrid.initialcondition!(soil, processes(soil), state)
end

end
