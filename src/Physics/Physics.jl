module Physics

using CryoGrid: Process, SubSurfaceProcess, CoupledProcesses, BoundaryProcess, Coupled, Layer, Top, Bottom, SubSurface
using CryoGrid: Event, DiscreteEvent, ContinuousEvent, ContinuousTrigger
using CryoGrid.Numerics
using CryoGrid.Utils

import CryoGrid

using Lazy: @>>
using Reexport
using Unitful

export volumetricfractions, waterice, liquidwater, partial

include("common.jl")
include("coupled.jl")
include("Boundaries/Boundaries.jl")
@reexport using .Boundaries
include("HeatConduction/HeatConduction.jl")
@reexport using .HeatConduction
include("Hydrology/Hydrology.jl")
@reexport using .Hydrology
include("Snow/Snow.jl")
@reexport using .Snow
include("Soils/Soils.jl")
@reexport using .Soils
include("SEB/SEB.jl")
@reexport using .SEB
include("Sources/Sources.jl")
@reexport using .Sources

end