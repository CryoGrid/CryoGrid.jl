module Diagnostics

using CryoGrid.Utils

using DimensionalData
using Statistics
using TimeSeries
using Unitful

export zero_annual_amplitude, permafrostdepth, active_layer_depth, thawdepth, mean_annual_ground_temperature

_check_array_dims(arr::AbstractDimArray) = @assert hasdim(arr, Z) && hasdim(arr, Ti) "array must have depth (Z) and time (Ti) dimensions"
"""
    zero_annual_amplitude(arr::AbstractDimArray{<:TempQuantity}; threshold=0.5u"°C")

Computes annual depth of zero amplitude (where `|max - min| < threshold`) and returns the
result for each year. Assumes `arr` to have dimensions `Ti` (time) and `Z` (depth) in any order.
"""
function zero_annual_amplitude(arr::AbstractDimArray{<:TempQuantity}; threshold=0.5u"K")
    _check_array_dims(arr)
    T = TimeArray(dims(arr, Ti).val, permutedims(arr, (Ti,Z)))
    annual_min = collapse(T, year, first, minimum)
    annual_max = collapse(T, year, first, maximum)
    z_inds = argmax(abs.(values(annual_max) .- values(annual_min)) .< threshold, dims=2)
    rebuild(arr, dims(arr, Z)[[i[2] for i in z_inds]], (Ti(timestamp(annual_max)),))
end
"""
    permafrostdepth(arr::AbstractDimArray{<:TempQuantity})

Computes depth of permafrost layer for all years, i.e. the closest depth to the surface at which
the maximum annual temperature is strictly `< 0°C`. Assumes `arr` to have dimensions `Ti` (time) and
`Z` (depth) in any order.
"""
function permafrostdepth(arr::AbstractDimArray{T}) where {T<:TempQuantity}
    _check_array_dims(arr)
    tarr = TimeArray(dims(arr, Ti).val, permutedims(arr, (Ti,Z)))
    annual_max = collapse(tarr, year, first, maximum)
    # find argmax (first one/"true" bit) of annual_max < 0
    z_inds = argmax(values(annual_max) .< 0u"°C", dims=2)
    rebuild(arr, dims(arr, Z)[[i[2] for i in z_inds]], (Ti(timestamp(annual_max)),))
end
"""
    active_layer_depth(arr::AbstractDimArray{<:TempQuantity})

Computes active layer depth annually. The active layer is assumed to extend from the
surface to the deepest point where the annual `max > 0` and `min < 0`, thus excluding any residual
thaw (non-cryotic) ground above the permafrost layer. Note that, in the unusual case where
even the surface temperature is always < 0 or > 0, the topmost (typically surface) depth will be returned.
Assumes `arr` to have dimensions `Ti` (time) and `Z` (depth) in any order.
"""
function active_layer_depth(arr::AbstractDimArray{<:TempQuantity})
    _check_array_dims(arr)
    # first get permafrost depth and select only grid cells above permafrost layer
    permadepth = maximum(permafrostdepth(arr))
    nonfrozen = arr[Z(Where(z -> z <= permadepth))]
    tarr = TimeArray(dims(arr, Ti).val, permutedims(nonfrozen, (Ti,Z)))
    # now get argmin, i.e. first depth where `(annual_min < 0) & (annual_max > 0)` is false
    annual_min = collapse(tarr, year, first, minimum)
    annual_max = collapse(tarr, year, first, maximum)
    z_inds = argmin((values(annual_min) .< 0u"°C") .& (values(annual_max) .> 0u"°C"), dims=2)
    # note that we subtract one to go up to the last cell where the condition held
    rebuild(arr, dims(arr, Z)[[max(1,i[2]-1) for i in z_inds]], (Ti(timestamp(annual_max)),))
end
"""
    thawdepth(arr::AbstractDimArray{<:TempQuantity})

Computes thaw depth (a.k.a freezing front) at all time steps. Assumes `arr` to have dimensions `Ti` (time)
and `Z` (depth) in any order.
"""
function thawdepth(arr::AbstractDimArray{<:TempQuantity})
    _check_array_dims(arr)
    T = permutedims(arr, (Z,Ti))
    z_inds = argmax((T .<= 0u"°C").data, dims=1)
    rebuild(arr, dims(arr, Z)[[i[1] for i in z_inds]], (dims(T,Ti),))
end
"""
    mean_annual_ground_temperature(arr::AbstractDimArray; upper_limit=0u"m", lower_limit=10u"m")

Computes mean annual ground temperature between `upper_limit` and `lower_limit`. Assumes `arr` to have dimensions
`Ti` (time) and `Z` (depth) in any order.
"""
function mean_annual_ground_temperature(arr::AbstractDimArray; upper_limit=0u"m", lower_limit=10u"m")
    _check_array_dims(arr)
    T = permutedims(ustrip.(arr), (Z,Ti))[Z(Between(upper_limit, lower_limit))]
    return mean(T, dims=1)[1,:]
end

end