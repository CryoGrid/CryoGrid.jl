module Physics

import CryoGrid: Process, CoupledProcesses, BoundaryProcess, Coupled, Layer, Top, Bottom, SubSurface
import CryoGrid: diagnosticstep!, initialcondition!, interact!, prognosticstep!, variables, callbacks, observe

using CryoGrid.Numerics
using CryoGrid.Utils

using Lazy: @>>
using Reexport
using Unitful

"""
    totalwater(sub::SubSurface, state)
    totalwater(sub::SubSurface, state)
    totalwater(sub::SubSurface, state, i)

Retrieves the total water content for the given layer at grid cell `i`, if specified.
Defaults to using the scalar total water content defined on layer `sub`.
"""
@inline totalwater(sub::SubSurface, state) = state.θw
@inline totalwater(sub::SubSurface, state, i) = Utils.getscalar(totalwater(sub, state), i)

include("coupled.jl")
include("Boundaries/Boundaries.jl")
include("HeatConduction/HeatConduction.jl")
include("WaterBalance/WaterBalance.jl")
include("Soils/Soils.jl")
include("SEB/SEB.jl")
include("Sources/Sources.jl")

@reexport using .Boundaries
@reexport using .HeatConduction
@reexport using .WaterBalance
@reexport using .Soils
@reexport using .SEB
@reexport using .Sources

end