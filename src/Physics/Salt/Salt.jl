module Salt

using CryoGrid
using CryoGrid.Heat
using CryoGrid.Heat: Temperature
using CryoGrid.Hydrology
using CryoGrid.Numerics
using CryoGrid.Soils
using CryoGrid.Utils

using FreezeCurves
using StaticArrays

export MarineSediment, SaltMassBalance, SaltProperties
include("salt_types.jl")

include("salt_diffusion.jl")

export SaltGradient
include("salt_bc.jl")

end
