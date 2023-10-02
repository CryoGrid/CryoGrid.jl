"""
    integrate(X::AbstractDimArray, grid::Grid{Edges}; upper_limit=0u"m", lower_limit=10u"m")

Integrates the quantity `X` over the given `grid`, which is assumed to be spatially alligned, i.e.
`length(grid) == length(dims(X,Z)) + 1` and `cells(grid) .≈ dims(X,Z)` are necessary preconditions.
"""
function integrate(X::AbstractDimArray, grid::Grid{Edges,G,<:DistQuantity}; upper_limit=0u"m", lower_limit=Inf*u"m") where G
    _check_arr_dims(X)
    X = permutedims(X, (Z,Ti))
    @assert length(grid) == size(X,1) + 1
    @assert all(cells(grid) .≈ dims(X,Z))
    Δx = Δ(grid)
    X_scaled = X.*Δx
    return sum(X_scaled[Z(Between(lower_limit, upper_limit))], dims=1)[1,:].*area(grid)
end

"""
    computejac(tile::Tile, u, p, t)

Helper function that computes the Jacobian of the given `tile` at `u` with parameters `p` and time `t`.
"""
function computejac(tile::Tile, u, p, t)
    J = Numerics.ForwardDiff.jacobian(u) do u
        du = zero(u)
        tile(du, u, p, t)
        return du
    end
    return J
end

"""
    build_dummy_state(grid::Grid, layer::Layer; t=0.0, with_units=true)

Collects variables defined on `layer` and initializes a dummy `state` named tuple
with all state variables initialized from the given `grid`. This is intended to be
used for unit tests and debugging in order to avoid the full-fledged construction
of a `Tile`/`Stratigraphy` and associated state types.
"""
function build_dummy_state(grid::Grid, layer::Layer; t=0.0, dt=1.0, with_units=true)
	vargrid(::OnGrid{Cells}, grid::Grid) = cells(grid)
	vargrid(::OnGrid{Edges}, grid::Grid) = edges(grid)
	maybeunits(var::Var) = with_units ? varunits(var) : Unitful.NoUnits
	vars = CryoGrid.variables(layer)
    # create delta vars
    dvars = map(CryoGrid.DVar, filter(isprognostic, vars))
    all_vars = tuple(vars..., dvars...)
	return (
		t = t,
        dt = dt*(with_units ? u"s" : 1),
        z = [first(grid)],
        Δz = [grid[end] - grid[1]],
		grid = with_units ? grid : Grid(ustrip.(grid)),
		map(v -> varname(v) => zeros(dimlength(vardims(v), length(grid)))*maybeunits(v), all_vars)...
	)
end
