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

export Soil, SoilParameterization, HomogeneousSoil
include("soil_types.jl")

export SoilProfile, porosity, mineral, organic
include("soil_methods.jl")

export SoilTexture
include("soil_texture.jl")

export MineralOrganic, soilcomponent
include("para/mineral_organic.jl")

export SURFEX
include("para/surfex.jl")

export RichardsEq
include("soil_water.jl")

export SoilThermalProperties
include("soil_heat.jl")

# water/heat coupling methods
include("soil_water_heat.jl")

end
