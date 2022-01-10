module WaterBalance

import CryoGrid: BoundaryStyle, diagnosticstep!, prognosticstep!, interact!, initialcondition!, variables

using ..Physics
using CryoGrid.Numerics
using CryoGrid.Numerics: nonlineardiffusion!

using IfElse

export SWRC, VanGenuchten

include("swrc.jl")

# TODO: implement water fluxes

end
