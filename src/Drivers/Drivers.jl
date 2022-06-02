"""
Module containing types and functions for enabling time-integration of a CryoGrid `Tile` model.
"""
module Drivers

using CryoGrid
using CryoGrid.Utils

using Reexport

export JacobianStyle, DefaultJac, TridiagJac, HeatOnlyTile

"""
    JacobianStyle

Trait for indicating Jacobian sparsity of a CryoGrid ODEProblem.
"""
abstract type JacobianStyle end
struct DefaultJac <: JacobianStyle end
struct TridiagJac <: JacobianStyle end
"""
    JacobianStyle(::Type{<:Tile})

Can be overriden/extended to specify Jacobian structure for specific `Tile`s.
"""
JacobianStyle(::Type{<:Tile}) = DefaultJac()
# Auto-detect Jacobian sparsity for problems with one or more heat-only layers.
# Note: This assumes that the processes/forcings on the boundary layers do not violate the tridiagonal structure!
# Unfortunately, the Stratigraphy type signature is a bit nasty to work with :(
const HeatOnlyTile = Tile{<:Stratigraphy{N,<:Tuple{TTop,Vararg{<:Union{<:StratComponent{<:SubSurface, <:CoupledProcesses{<:Tuple{<:Heat}}},TBot}}}}} where {N,TTop,TBot}
JacobianStyle(::Type{<:HeatOnlyTile}) = TridiagJac()

# DiffEq/SciML driver (possibly should be a soft dependency with Requires.jl)
include("DiffEq/DiffEq.jl")
@reexport using .DiffEq

end
