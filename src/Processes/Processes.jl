module Processes

import CryoGrid.Interface: diagnosticstep!, initialcondition!, interact!, prognosticstep!, variables

using CryoGrid.Interface
using CryoGrid.Layers
using CryoGrid.Numerics
using CryoGrid.Utils

using Lazy: @>>
using Reexport
using Unitful

include("systems.jl")
include("Boundaries/Boundaries.jl")
include("Water/Water.jl")
include("HeatConduction/HeatConduction.jl")
include("SEB/SEB.jl")
include("Sources/Sources.jl")

@reexport using .Boundaries
@reexport using .HeatConduction
@reexport using .SEB
@reexport using .Sources

end