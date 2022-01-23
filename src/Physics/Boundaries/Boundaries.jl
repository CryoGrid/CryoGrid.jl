module Boundaries

import CryoGrid: BoundaryProcess, BoundaryStyle, Dirichlet, Neumann, Top
import CryoGrid: variables, boundaryvalue

using CryoGrid.Numerics
using CryoGrid.Utils

using Base: @propagate_inbounds
using ConstructionBase
using Dates
using Flatten
using Interpolations
using ModelParameters
using Parameters
using TimeSeries
using Unitful

import Flatten: flattenable

export BoundaryEffect, Damping
include("effects.jl")
export Forcing, TimeSeriesForcing, ForcingData
include("forcing.jl")
export ConstantBC, PeriodicBC, Bias
include("bc.jl")

end
