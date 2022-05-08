module Boundaries

import CryoGrid: SubSurfaceProcess, BoundaryProcess, BoundaryStyle, Dirichlet, Neumann, Top
import CryoGrid: variables, boundaryvalue

using CryoGrid.Numerics
using CryoGrid.Utils

using ConstructionBase
using Dates
using Unitful

import Flatten: flattenable

"""
    BoundaryEffect

Base type for boundary "effects" which modify boundary conditions based on some
given parameterization.
"""
abstract type BoundaryEffect end

export BoundaryEffect

export ConstantBC, PeriodicBC, Bias
include("bc.jl")

end
