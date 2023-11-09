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

export Snowpack, SnowMassBalance, SnowBC
include("snow_types.jl")

# dynamic snow mass
include("snow_mass.jl")

include("snow_heat.jl")

# single-layer "bulk" snow scheme
include("snow_bulk.jl")

export swe, snowdensity, snowdepth
include("snow_methods.jl")

include("snow_bc.jl")

end