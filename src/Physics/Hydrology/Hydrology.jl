module Hydrology

import CryoGrid: BoundaryStyle, diagnosticstep!, prognosticstep!, interact!, initialcondition!, variables

using ..Physics
using CryoGrid.Numerics
using CryoGrid.Numerics: nonlineardiffusion!

using IfElse
using ModelParameters
using Unitful

include("water_bucket.jl")

end
