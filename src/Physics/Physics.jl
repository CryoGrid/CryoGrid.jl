module Physics

import CryoGrid: Process, CompoundProcess, Coupled, Layer, Top, Bottom, SubSurface
import CryoGrid: diagnosticstep!, initialcondition!, interact!, prognosticstep!, variables, observe

using CryoGrid.Layers
using CryoGrid.Numerics
using CryoGrid.Utils

using Lazy: @>>
using Reexport
using Unitful

include("compound.jl")
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