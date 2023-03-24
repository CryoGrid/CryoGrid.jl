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

export SnowMassBalance
include("snow_mass.jl")

export Snowpack, SnowProperties
export swe, snowdensity
include("snowpack.jl")

export Snowfall
include("snow_bc.jl")

include("snow_bulk.jl")

end