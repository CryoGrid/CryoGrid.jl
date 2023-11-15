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

# Local alias for Heat EnthalpyBased type
const EnthalpyBased = Heat.EnthalpyBased

include("snow_types.jl")

export Snowpack, SnowMassBalance, SnowBC
include("snowpack.jl")

include("snow_mass.jl")

include("snow_heat.jl")

# single-layer "bulk" snow scheme
include("snow_bulk.jl")

export swe, snowdensity, snowdepth, snowfall
export accumulation!, ablation!
include("snow_methods.jl")

include("snow_bc.jl")

end