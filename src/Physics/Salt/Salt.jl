module Salt

using CryoGrid
using CryoGrid.Heat
using CryoGrid.Heat: TemperatureBased
using CryoGrid.Hydrology
using CryoGrid.Numerics
using CryoGrid.Soils
using CryoGrid.Utils

using ForwardDiff
using FreezeCurves
using StaticArrays

export SalineGround, SaltMassBalance, SaltProperties
include("salt_types.jl")

include("salt_diffusion.jl")

export SaltBC, SaltGradient, SaltHeatBC
include("salt_bc.jl")

export SedimentCompactionInitializer
include("salt_init.jl")

end
