module LiteImplicit

using CryoGrid
using CryoGrid.Numerics
using CryoGrid.Utils

using DataStructures
using DiffEqBase, DiffEqCallbacks
using Interpolations

import CommonSolve
import SciMLBase

export LiteImplicitEuler
include("cglite_types.jl")

const CGLiteIntegrator = CryoGridIntegrator{T} where {T<:LiteImplicitEuler}

include("cglite_step.jl")

end
