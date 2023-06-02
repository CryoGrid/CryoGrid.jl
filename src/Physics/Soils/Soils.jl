module Soils

using CryoGrid
using CryoGrid.Heat
using CryoGrid.Hydrology
using CryoGrid.Numerics
using CryoGrid.Utils

using IfElse
using IntervalSets
using ForwardDiff
using ModelParameters
using Reexport: @reexport
using Setfield
using StaticArrays
using Unitful
using UnPack

# Re-export FreezeCurves package
@reexport using FreezeCurves
@reexport using FreezeCurves.Solvers

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
