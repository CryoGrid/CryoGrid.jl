module Analysis

using CryoGrid.Utils

using DimensionalData
using TimeSeries

_check_array_dims(arr::AbstractDimArray) = @assert hasdim(arr, Z) && hasdim(arr, Ti) "array must have depth (Z) and time (Ti) dimensions"

"""
    zero_annual_amplitude(arr::AbstractDimArray{<:TempQuantity}; threshold=0.5u"°C")

Computes annual depth of zero amplitude (where `|max - min| < threshold`) and returns the
result for each year. Assumes `arr` to have dimensions `Ti` (time) and `Z` (depth) in any order.
"""
function zero_annual_amplitude(arr::AbstractDimArray{<:TempQuantity}; threshold=0.5u"°C")
    _check_array_dims(arr)
    tarr = TimeArray(dims(arr, Ti).val, permutedims(arr, (Ti,Z)))
    annual_min = collapse(tarr, year, last, minimum)
    annual_max = collapse(tarr, year, last, maximum)
    z_inds = argmax(abs.(values(annual_min) .- values(annual_max)) .< threshold, dims=2)
    rebuild(arr, dims(x, Z)[[i[2] for i in z_inds]], (Ti(timestamp(annual_max)),))
end

end