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

export Snowpack, SnowMassBalance, DynamicSnowMassBalance, PrescribedSnowMassBalance, SnowBC
include("snow_types.jl")

export swe, snowdensity
include("snow_methods.jl")

# snow ablation/accumulation schemes
include("snow_mass.jl")

# snow density schemes
include("snow_density.jl")

# bulk snow scheme
include("snow_bulk.jl")

end