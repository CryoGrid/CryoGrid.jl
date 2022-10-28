module Physics

import CryoGrid

using CryoGrid
using CryoGrid.Numerics
using CryoGrid.Utils

using IfElse
using Reexport
using Unitful

export volumetricfractions, partial

include("common.jl")
include("steplimiters.jl")
include("Boundaries/Boundaries.jl")
include("Hydrology/Hydrology.jl")
include("Heat/Heat.jl")
include("Snow/Snow.jl")
include("Soils/Soils.jl")
include("SEB/SEB.jl")
include("Sources/Sources.jl")

@reexport using .Boundaries
@reexport using .Heat
@reexport using .Hydrology
@reexport using .Snow
@reexport using .Soils
@reexport using .SEB
@reexport using .Sources

end