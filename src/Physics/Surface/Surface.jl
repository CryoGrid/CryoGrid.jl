module Surface

using CryoGrid
using CryoGrid.Heat
using CryoGrid.Hydrology
using CryoGrid.Snow
using CryoGrid.Soils
using CryoGrid.Numerics
using CryoGrid.Utils

using NonlinearSolve
using Setfield
using StaticArrays: @SVector
using Unitful

import Flatten

# Surface energy balance
export SurfaceEnergyBalance, SEBParams
include("SEB/seb.jl")

export SurfaceWaterBalance
include("SWB/swb.jl")

end
