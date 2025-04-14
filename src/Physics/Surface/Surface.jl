module Surface

using CryoGrid
using CryoGrid.Heat
using CryoGrid.Hydrology
using CryoGrid.Snow
using CryoGrid.Soils
using CryoGrid.Numerics
using CryoGrid.Utils

using Setfield
using NonlinearSolve
using StaticArrays: @SVector
using Unitful

import Flatten

# Surface energy balance
export SurfaceEnergyBalance, SEBParams
export relative_to_specific_humidity, vapor_pressure
include("SEB/seb.jl")

export SurfaceWaterBalance
include("SWB/swb.jl")

end
