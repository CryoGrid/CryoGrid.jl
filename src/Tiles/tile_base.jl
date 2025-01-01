"""
    AbstractTile{iip}

Base type for 1D tiles. `iip` is a boolean value that indicates, if true,
whether the model operates on state variables in-place (overwriting arrays) or
if false, out-of-place (copying arrays). Current only in-place is supported.
"""
abstract type AbstractTile{iip} end

"""
    (tile::AbstractTile{true})(du,u,p,t,dt=1.0)
    (tile::AbstractTile{false})(u,p,t,dt=1.0)

Invokes `computeprognostic!` on `tile` to compute the time derivative du/dt.
"""
(tile::AbstractTile{true})(du, u, p, t, dt=1.0) = computeprognostic!(tile,du,u,p,t,dt)
(tile::AbstractTile{false})(u, p, t, dt=1.0) = computefluxes(tile,u,p,t,dt)

"""
    isinplace(tile::AbstractTile{iip}) where {iip}

Returns true if `tile` uses in-place mode, false if out-of-place.
"""
isinplace(::AbstractTile{iip}) where {iip} = iip

"""
    computeprognostic!(::T,du,u,p,t) where {T<:AbstractTile}

In-place update function for tile `T`. Computes du/dt and stores the result in `du`.
"""
computeprognostic!(::T, du, u, p, t, dt=1.0) where {T<:AbstractTile} = error("no implementation of in-place computeprognostic! for $T")

"""
    computefluxes(::T,u,p,t) where {T<:AbstractTile}

Out-of-place update function for tile `T`. Computes and returns du/dt as vector with same size as `u`.
"""
computefluxes(::T, u, p, t, dt=1.0) where {T<:AbstractTile} = error("no implementation of out-of-place evaluate for $T")
