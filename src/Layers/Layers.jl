module Layers

import CryoGrid.Interface: SubSurface, initialcondition!, variables

using CryoGrid.Interface
using CryoGrid.Numerics
using CryoGrid.Utils

using DimensionalData
using IntervalSets
using Parameters
using Unitful

include("soil.jl")

end