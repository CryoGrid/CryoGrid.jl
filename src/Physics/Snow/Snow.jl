module Snow

using CryoGrid
using CryoGrid: ContinuousEvent, Increasing, Decreasing # for events/callbacks
using CryoGrid.Physics
using CryoGrid.Numerics
using CryoGrid.Utils

import CryoGrid
import CryoGrid.InputOutput
import CryoGrid.Physics
import CryoGrid.Physics.Heat

using IfElse
using ModelParameters
using Unitful
using UnPack

SnowThermalProperties = Heat.ThermalProperties

export Snowpack, SnowProperties, SnowMassBalance
include("types.jl")

include("methods.jl")

export Snowfall
include("snow_bc.jl")

include("snow_bulk.jl")

end