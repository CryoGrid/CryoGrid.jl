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

import ConstructionBase

# aliases for heat formulations in Heat module
const TemperatureBased = Heat.TemperatureBased
const EnthalpyBased = Heat.EnthalpyBased
const EnthalpyImplicit = Heat.EnthalpyImplicit

export Ground, AbstractGround
include("ground.jl")

export Soil, SoilParameterization, Heterogeneous
include("soil_types.jl")

export SoilProfile, porosity, mineral, organic
include("soil_methods.jl")

export SoilTexture
include("soil_texture.jl")

export Heterogeneous, SimpleSoil, SURFEXSoil
export soilcomponent
include("para/soil_para.jl")

export RichardsEq
include("soil_water.jl")

export SoilThermalProperties
include("soil_heat.jl")

# water/heat coupling methods
include("soil_water_heat.jl")

end
