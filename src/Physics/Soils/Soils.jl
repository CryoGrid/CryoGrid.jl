module Soils

using CryoGrid
using CryoGrid.Numerics
using CryoGrid.Utils

using IfElse
using IntervalSets
using ForwardDiff
using FreezeCurves
using FreezeCurves.Solvers
using ModelParameters
using Setfield
using StaticArrays
using Unitful
using UnPack

import CryoGrid
import CryoGrid.InputOutput
import CryoGrid.Heat
import CryoGrid.Hydrology

# from FreezeCurves.jl
export SFCC, PainterKarra, DallAmico, DallAmicoSalt, Westermann, McKenzie
export VanGenuchten, BrooksCorey

# aliases for heat formulations in Heat module
const Temperature = Heat.Temperature
const Enthalpy = Heat.Enthalpy
const EnthalpyImplicit = Heat.EnthalpyImplicit

export Soil, HomogeneousSoil, SoilParameterization, SoilProfile, SoilProperties
export soilproperties, porosity, mineral, organic
include("soil.jl")

export HomogeneousMixture, MineralSediment, soilcomponent
include("soil_para.jl")

export RichardsEq
include("soil_water.jl")

export SoilThermalProperties
include("soil_heat.jl")

# water/heat coupling methods
include("soil_water_heat.jl")

end
