module Snow

using CryoGrid: Top, SubSurfaceProcess, BoundaryProcess, BoundaryStyle, Dirichlet
using CryoGrid.InputOutput: Forcing
using CryoGrid.Physics.HeatConduction
using CryoGrid.Numerics
using CryoGrid.Utils

import CryoGrid
import CryoGrid.Physics
import CryoGrid.Physics.HeatConduction

using Unitful

struct SnowFree <: CryoGrid.Callback end
CryoGrid.CallbackStyle(::Type{SnowFree}) = CryoGrid.Continuous()

export Snowpack
include("snowpack.jl")
export SnowMassBalance, Snowfall
include("processes/snow_processes.jl")

end