module Physics

import CryoGrid: Process, CoupledProcesses, BoundaryProcess, Coupled, Layer, Top, Bottom, SubSurface
import CryoGrid: diagnosticstep!, initialcondition!, interact!, prognosticstep!, variables, callbacks, observe

using CryoGrid.Numerics
using CryoGrid.Utils

using Lazy: @>>
using Reexport
using Unitful

export waterice

include("common.jl")
include("coupled.jl")
include("Boundaries/Boundaries.jl")
include("HeatConduction/HeatConduction.jl")
include("WaterBalance/WaterBalance.jl")
include("Snow/Snow.jl")
include("Soils/Soils.jl")
include("SEB/SEB.jl")
include("Sources/Sources.jl")

@reexport using .Boundaries
@reexport using .HeatConduction
@reexport using .WaterBalance
@reexport using .Snow
@reexport using .Soils
@reexport using .SEB
@reexport using .Sources

end