module Processes

import CryoGrid.Interface: diagnosticstep!, initialcondition!, interact!, prognosticstep!, variables

using CryoGrid.Interface
using CryoGrid.Layers
using CryoGrid.Numerics
using CryoGrid.Utils

using Lazy: @>>
using Reexport
using Unitful

include("boundaries.jl")
include("systems.jl")
include("HeatConduction/HeatConduction.jl")
include("SEB/SEB.jl")

@reexport using .HeatConduction
@reexport using .SEB
@reexport using .Systems

end