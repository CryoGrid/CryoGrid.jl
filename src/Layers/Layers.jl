module Layers

import CryoGrid: SubSurface
import CryoGrid: initialcondition!, variables

using CryoGrid.Numerics
using CryoGrid.Utils

using DimensionalData
using IntervalSets
using Parameters
using Unitful

export Soil, SoilProperties, SoilProfile, SoilParams, SoilType, Sand, Silt, Clay, SoilParameterization, BySoilProperties

include("soil.jl")

end