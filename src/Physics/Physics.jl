module Physics

import CryoGrid: Process, CoupledProcesses, BoundaryProcess, Coupled, Layer, Top, Bottom, SubSurface
import CryoGrid: diagnosticstep!, initialcondition!, interact!, prognosticstep!, variables, callbacks, observe

using CryoGrid.Numerics
using CryoGrid.Utils

using Lazy: @>>
using Reexport
using Unitful

"""
    waterice(sub::SubSurface, state)
    waterice(sub::SubSurface, state, i)

Retrieves the total water content (water + ice) for the given layer at grid cell `i`, if specified.
Defaults to retrieving the state variable `θwi`. 
"""
@inline waterice(::SubSurface, state) = state.θwi
@inline waterice(sub::SubSurface, state, i) = Utils.getscalar(waterice(sub, state), i)

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