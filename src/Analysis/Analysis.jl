module Analysis

using CryoGrid.Utils

using DimensionalData
using TimeSeries
using Unitful

_check_array_dims(arr::AbstractDimArray) = @assert hasdim(arr, Z) && hasdim(arr, Ti) "array must have depth (Z) and time (Ti) dimensions"

"""
    zero_annual_amplitude(arr::AbstractDimArray{<:TempQuantity}; threshold=0.5u"°C")

Computes annual depth of zero amplitude (where `|max - min| < threshold`) and returns the
result for each year. Assumes `arr` to have dimensions `Ti` (time) and `Z` (depth) in any order.
"""
function zero_annual_amplitude(arr::AbstractDimArray{<:TempQuantity}; threshold=0.5u"°C")
    _check_array_dims(arr)
    tarr = TimeArray(dims(arr, Ti).val, permutedims(arr, (Z,Ti)))
    annual_min = collapse(tarr, year, last, minimum)
    annual_max = collapse(tarr, year, last, maximum)
    z_inds = argmax(abs.(values(annual_min) .- values(annual_max)) .< threshold, dims=1)
    rebuild(arr, dims(x, Z)[[i[1] for i in z_inds]], (Ti(timestamp(annual_max)),))
end

"""
    permafrostdepth(arr::AbstractDimArray{<:TempQuantity})

Computes depth of permafrost layer for all years, i.e. the closest depth to the surface at which
the maximum annual temperature is <= 0°C. Assumes `arr` to have dimensions `Ti` (time) and
`Z` (depth) in any order.
"""
function permafrostdepth(arr::AbstractDimArray{<:TempQuantity})
    _check_array_dims(arr)
    tarr = TimeArray(dims(arr, Ti).val, permutedims(arr, (Z,Ti)))
    annual_max = collapse(tarr, year, last, maximum)
    # find argmax (first one/"true" bit) of annual_max <= 0
    z_inds = argmax(values(annual_max .<= 0), dims=1)
    rebuild(arr, dims(x, Z)[[i[1] for i in z_inds]], (Ti(timestamp(annual_max)),))
end

"""
    active_layer_thickness(arr::AbstractDimArray{<:TempQuantity})

Computes active layer thickness for all years. The active layer is assumed to extend from the
surface to the deepest point where the annual `max > 0` and `min < 0`, thus excluding any residual
thaw (non-cryotic) ground above the permafrost layer. Note that, in the unusual case where
even the surface temperature is always < 0 or > 0, the topmost (typically surface) depth will be returned.
Assumes `arr` to have dimensions `Ti` (time) and `Z` (depth) in any order.
"""
function active_layer_thickness(arr::AbstractDimArray{<:TempQuantity})
    _check_array_dims(arr)
    # first get permafrost depth and select only grid cells above permafrost layer
    perma_depth = permafrostdepth(arr).data
    nonfrozen = arr[Z(perma_depth)]
    tarr = TimeArray(dims(arr, Ti).val, permutedims(nonfrozen, (Z,Ti)))
    # now get argmin, i.e. first depth where `(annual_min < 0) & (annual_max > 0)` is false
    annual_min = collapse(tarr, year, last, minimum)
    annual_max = collapse(tarr, year, last, maximum)
    z_inds = argmin((values(annual_min) .< 0) .& (values(annual_max) .> 0), dims=1)
    # note that we subtract one to go up to the last cell where the condition held
    rebuild(arr, dims(x, Z)[[min(1,i[1]-1) for i in z_inds]], (Ti(timestamp(annual_max)),))
end

"""
    thawdepth(arr::AbstractDimArray{<:TempQuantity})

Computes thaw depth (a.k.a freezing front) at all time steps. Assumes `arr` to have dimensions `Ti` (time)
and `Z` (depth) in any order.
"""
function thawdepth(arr::AbstractDimArray{<:TempQuantity})
    _check_array_dims(arr)
    T = permutedims(arr, (Z,Ti))
    z_inds = argmax(values(T .<= 0), dims=1)
    rebuild(arr, dims(x, Z)[[i[1] for i in z_inds]], (Ti(timestamp(T)),))
end

"""
    mean_annual_ground_temperature(arr::AbstractDimArray{<:TempQuantity}; upper_limit=0u"m", lower_limit=10u"m")

Computes mean annual ground temperature between `upper_limit` and `lower_limit`. Assumes `arr` to have dimensions
`Ti` (time) and `Z` (depth) in any order.
"""
function mean_annual_ground_temperature(arr::AbstractDimArray{<:TempQuantity}; upper_limit=0u"m", lower_limit=10u"m")
    _check_array_dims(arr)
    T = permutedims(arr, (Ti,Z))[Z(Between(upper_limit, lower_limit))]
    return pemutedims(mean(T, dims=2), typeof.(dims(arr)))
end

end