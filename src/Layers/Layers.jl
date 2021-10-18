module Layers

import CryoGrid: SubSurface
import CryoGrid: initialcondition!, variables

using CryoGrid.Numerics
using CryoGrid.Utils

using DimensionalData
using IntervalSets
using ModelParameters
using Parameters
using Unitful

export Soil, SoilComposition, SoilProfile, SoilType, Sand, Silt, Clay

include("soil.jl")

end