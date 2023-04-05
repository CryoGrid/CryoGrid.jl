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
        du = similar(u)
        Strat.step!(tile, du, u, p, t)
        return du
    end
    return J
end
