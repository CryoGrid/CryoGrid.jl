"""
Module containing types and functions for enabling time-integration of a CryoGrid land model.
"""
module Drivers

using CryoGrid: Land, SubSurface, CompoundProcess
using CryoGrid.Numerics: varname
using CryoGrid.Physics: Heat
using CryoGrid.Utils: convert_tspan

using Dates
using ModelParameters
using LinearAlgebra
using Reexport

import CryoGrid.Land: LandModel, Stratigraphy, StratComponent, componenttypes, boundaries, withaxes, getstate, getstates, getvar

export DefaultJac, TridiagJac, JacobianStyle, HeatOnlyLandModel

"""
    JacobianStyle

Trait for indicating Jacobian sparsity of a CryoGrid ODEProblem.
"""
abstract type JacobianStyle end
struct DefaultJac <: JacobianStyle end
struct TridiagJac <: JacobianStyle end
"""
    JacobianStyle(::Type{<:LandModel})

Can be overriden/extended to specify Jacobian structure for specific `LandModel`s.
"""
JacobianStyle(::Type{<:LandModel}) = DefaultJac()
# Auto-detect Jacobian sparsity for problems with one or more heat-only layers.
# Note: This assumes that the processes/forcings on the boundary layers do not violate the tridiagonal structure!
# Unfortunately, the Stratigraphy type signature is a bit nasty to work with :(
const HeatOnlyLandModel = LandModel{<:Stratigraphy{N,<:Tuple{TTop,Vararg{<:Union{<:StratComponent{<:SubSurface, <:CompoundProcess{<:Tuple{<:Heat}}},TBot}}}}} where {N,TTop,TBot}
JacobianStyle(::Type{<:HeatOnlyLandModel}) = TridiagJac()

# DiffEq/SciML driver
export CryoGridProblem
include("diffeq.jl")

end
