module Snow

using CryoGrid
using CryoGrid: ContinuousEvent, Increasing, Decreasing # for events/callbacks
using CryoGrid.Numerics
using CryoGrid.Utils

import CryoGrid
import CryoGrid.InputOutput
import CryoGrid.Heat

using IfElse
using ModelParameters
using Unitful
using UnPack

export SnowMassBalance, SnowBC
include("snow_mass.jl")

export Snowpack, SnowProperties
export swe, snowdensity
include("snowpack.jl")

include("snow_bulk.jl")

end