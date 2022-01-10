module Diagnostics

using CryoGrid.Numerics
using CryoGrid.Strat
using CryoGrid.Utils

using DimensionalData
using Statistics
using TimeSeries
using Unitful

export zero_annual_amplitude, permafrosttable, permafrostbase, thawdepth, active_layer_thickness, mean_annual_ground_temperature

include("spinup.jl")

export spinup

_check_arr_dims(T::AbstractDimArray) = @assert hasdim(T, Z) && hasdim(T, Ti) "Tay must have depth (Z) and time (Ti) dimensions"
"""
    zero_annual_amplitude(T::AbstractDimArray{<:TempQuantity}; threshold=0.5u"°C")

Computes annual depth of zero amplitude (where `|max - min| < threshold`) and returns the
result for each year. Assumes `T` to have dimensions `Ti` (time) and `Z` (depth) in any order.
"""
function zero_annual_amplitude(T::AbstractDimArray{<:TempQuantity}; threshold=0.5u"K")
    _check_arr_dims(T)
    Tarr = TimeArray(dims(T, Ti).val, permutedims(T, (Ti,Z)))
    annual_min = collapse(Tarr, year, first, minimum)
    annual_max = collapse(Tarr, year, first, maximum)
    z_inds = vec(argmax(abs.(values(annual_max) .- values(annual_min)) .< threshold, dims=2))
    rebuild(T, collect(dims(T, Z)[[i[2] for i in z_inds]]), (Ti(timestamp(annual_max)),))
end
"""
    permafrosttable(T::AbstractDimArray{<:TempQuantity})

Computes depth of permafrost table for all years, i.e. the closest depth to the surface at which
the maximum annual temperature is strictly `< 0°C`. Assumes `T` to have dimensions `Ti` (time) and
`Z` (depth) in any order.
"""
function permafrosttable(T::AbstractDimArray{<:TempQuantity})
    _check_arr_dims(T)
    Tarr = TimeArray(dims(T, Ti).val, permutedims(T, (Ti,Z)))
    annual_max = collapse(Tarr, year, first, maximum)
    # find argmax (first one/"true" bit) of annual_max < 0
    z_inds = vec(argmax(values(annual_max) .< 0u"°C", dims=2))
    rebuild(T, collect(dims(T, Z)[[i[2] for i in z_inds]]), (Ti(timestamp(annual_max)),))
end
"""
    permafrostbase(T::AbstractDimArray{<:TempQuantity})

Computes depth of permafrost base for all years, i.e. the closest depth to the "bottom" at which
the maximum annual temperature is strictly `< 0°C`. Assumes `T` to have dimensions `Ti` (time) and
`Z` (depth) in any order.
"""
function permafrostbase(T::AbstractDimArray{<:TempQuantity})
    _check_arr_dims(T)
    Tarr = TimeArray(dims(T, Ti).val, permutedims(T, (Ti,Z)))
    annual_max = collapse(Tarr, year, first, maximum)
    # find argmax (first one/"true" bit) of annual_max < 0
    z_inds = vec(size(Tarr,2) .- [i[2] for i in argmax(reverse(values(annual_max) .< 0u"°C", dims=2), dims=2)])
    rebuild(T, collect(dims(T, Z)[[i for i in z_inds]]), (Ti(timestamp(annual_max)),))
end
"""
    thawdepth(T::AbstractDimArray{<:TempQuantity})

Computes thaw depth (a.k.a freezing front) at all time steps. Assumes `T` to have dimensions `Ti` (time)
and `Z` (depth) in any order.
"""
function thawdepth(T::AbstractDimArray{<:TempQuantity})
    _check_arr_dims(T)
    T = permutedims(T, (Z,Ti))
    z_inds = vec(argmax((T .<= 0u"°C").data, dims=1))
    rebuild(T, collect(dims(T, Z)[[i[1] for i in z_inds]]), (dims(T,Ti),))
end
"""
    active_layer_thickness(T::AbstractDimArray{<:TempQuantity})

Computes active layer thickness annually. The active layer thickness is defined here as the maximum thaw
depth throughout the calendar year. Assumes `T` to have dimensions `Ti` (time) and `Z` (depth) in any order.
"""
function active_layer_thickness(T::AbstractDimArray{<:TempQuantity})
    _check_arr_dims(T)
    thaw_depth = thawdepth(T)
    Tarr = TimeArray(dims(T, Ti).val, thaw_depth)
    annual_max = collapse(Tarr, year, first, maximum)
    rebuild(T, values(annual_max), (Ti(timestamp(annual_max)),))
end
"""
    mean_annual_ground_temperature(T::AbstractDimArray; upper_limit=0u"m", lower_limit=10u"m")

Computes mean annual ground temperature between `upper_limit` and `lower_limit`. Assumes `T` to have dimensions
`Ti` (time) and `Z` (depth) in any order.
"""
function mean_annual_ground_temperature(T::AbstractDimArray{<:TempQuantity}; upper_limit=0u"m", lower_limit=10u"m")
    _check_arr_dims(T)
    T = permutedims(ustrip.(T), (Z,Ti))[Z(Between(lower_limit, upper_limit))]
    return mean(T, dims=1)[1,:]
end
"""
    integrate(X::AbstractDimArray, grid::Grid{Edges}; upper_limit=0u"m", lower_limit=10u"m")

Integrates the quantity `X` over the given `grid`, which is assumed to be spatially alligned, i.e.
`length(grid) == length(dims(X,Z)) + 1` and `cells(grid) .≈ dims(X,Z)` are necessary preconditions.
"""
function integrate(X::AbstractDimArray, grid::Grid{Edges,G,<:DistQuantity}; upper_limit=0u"m", lower_limit=10u"m") where G
    _check_arr_dims(X)
    X = permutedims(X, (Z,Ti))
    @assert length(grid) == size(X,1) + 1
    @assert all(cells(grid) .≈ dims(X,Z))
    Δx = Δ(grid)
    X_scaled = X.*Δx
    return sum(X_scaled[Z(Between(lower_limit, upper_limit))], dims=1)[1,:].*area(grid)
end

end