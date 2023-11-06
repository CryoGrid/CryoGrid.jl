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

export Snowpack, SnowMassBalance, DynamicSnowMassBalance, PrescribedSnowMassBalance
include("snow_types.jl")

export SnowBC
include("snow_bc.jl")

# prescribed snow mass
include("snow_mass_prescribed.jl")

# dynamic snow mass
include("snow_mass_dynamic.jl")

# snow density schemes
include("snow_density.jl")

# single-layer "bulk" snow scheme
include("bulk/snow_bulk.jl")

export swe, snowdensity, snowdepth
include("snow_methods.jl")

end