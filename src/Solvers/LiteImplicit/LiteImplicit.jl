module LiteImplicit

using CryoGrid
using CryoGrid.Numerics
using CryoGrid.Utils

using DiffEqBase, DiffEqCallbacks
using Interpolations

import SciMLBase

export LiteImplicitEuler
include("cglite_types.jl")

const CGLiteIntegrator = CryoGridIntegrator{T} where {T<:LiteImplicitEuler}

include("step.jl")

end
