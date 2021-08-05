module Water

import CryoGrid.Interface: BoundaryStyle, diagnosticstep!, prognosticstep!, interact!, initialcondition!, variables

using ..Processes
using CryoGrid.Interface
using CryoGrid.Numerics
using CryoGrid.Numerics: nonlineardiffusion!

using IfElse

export SWRC, VanGenuchten

include("swrc.jl")

# TODO: implement water fluxes

end
