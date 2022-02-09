module Boundaries

import CryoGrid: BoundaryProcess, BoundaryStyle, Dirichlet, Neumann, Top
import CryoGrid: variables, boundaryvalue

using CryoGrid.Numerics
using CryoGrid.Utils

using ConstructionBase
using Dates
using Parameters
using Unitful

import Flatten: flattenable

export BoundaryEffect, Damping
include("effects.jl")
export ConstantBC, PeriodicBC, Bias
include("bc.jl")

end
