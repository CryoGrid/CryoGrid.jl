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

export Soil, SoilParameterization
include("types.jl")

export SoilProfile, SoilProperties, soilproperties, porosity, mineral, organic
include("methods.jl")

export HomogeneousMixture, MineralSediment, soilcomponent
include("soil_para.jl")

export RichardsEq
include("soil_water.jl")

export SoilThermalProperties
include("soil_heat.jl")

include("soil_water_heat_coupled.jl")

end
