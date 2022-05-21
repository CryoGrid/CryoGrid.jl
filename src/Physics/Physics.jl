module Physics

using CryoGrid: Process, SubSurfaceProcess, CoupledProcesses, BoundaryProcess, Coupled, Layer, Top, Bottom, SubSurface
using CryoGrid: Event, DiscreteEvent, ContinuousEvent, ContinuousTrigger
import CryoGrid: diagnosticstep!, initialcondition!, interact!, prognosticstep!, variables, events, observe

using CryoGrid.Numerics
using CryoGrid.Utils

using Lazy: @>>
using Reexport
using Unitful

export volumetricfractions, waterice

include("common.jl")
include("coupled.jl")
include("Boundaries/Boundaries.jl")
include("HeatConduction/HeatConduction.jl")
include("Hydrology/Hydrology.jl")
include("Snow/Snow.jl")
include("Soils/Soils.jl")
include("SEB/SEB.jl")
include("Sources/Sources.jl")

@reexport using .Boundaries
@reexport using .HeatConduction
@reexport using .Hydrology
@reexport using .Snow
@reexport using .Soils
@reexport using .SEB
@reexport using .Sources

end