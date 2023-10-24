"""
    JacobianStyle

Trait for indicating Jacobian sparsity of a CryoGrid ODEProblem.
"""
abstract type JacobianStyle end
struct DefaultJac <: JacobianStyle end
struct TridiagJac <: JacobianStyle end
"""
    JacobianStyle(::AbstractTile)

Can be overriden/extended to specify Jacobian structure for specific `Tile`s.
"""
function JacobianStyle(tile::AbstractTile)
    prognostic_grid_vars = filter(isongrid, filter(isprognostic, CryoGrid.variables(tile)))
    if length(prognostic_grid_vars) == 1
        # Auto-detect Jacobian sparsity for problems with one on-grid prognostic variable
        return TridiagJac()
    else
        return DefaultJac()
    end
end
