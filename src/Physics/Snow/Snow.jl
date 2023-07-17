module Snow

using CryoGrid
using CryoGrid: ContinuousEvent, Increasing, Decreasing # for events/callbacks
using CryoGrid.InputOutput
using CryoGrid.Heat
using CryoGrid.Hydrology
using CryoGrid.Numerics
using CryoGrid.Utils

using IfElse
using ModelParameters
using Unitful
using UnPack

export Snowpack, SnowMassBalance, SnowProperties, SnowMassProperties, SnowBC
include("snowpack.jl")

export swe, snowdensity
include("snow_methods.jl")

# bulk snow scheme
include("snow_bulk.jl")

end