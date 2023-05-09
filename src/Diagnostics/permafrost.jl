"""
    zero_annual_amplitude(T::AbstractDimArray{<:TempQuantity}; threshold=0.5u"K")

Computes annual depth of zero amplitude (where `|max - min| < threshold`) and returns the
result for each year. Assumes `T` to have dimensions `Ti` (time) and `Z` (depth) in any order.
"""
function zero_annual_amplitude(T::AbstractDimArray{<:TempQuantity}; threshold=0.5u"K")
    _check_arr_dims(T)
    Tarr = TimeArray(dims(T, Ti).val, permutedims(T, (Ti,Z)))
    annual_min = collapse(Tarr, year, first, minimum)
    annual_max = collapse(Tarr, year, first, maximum)
    T_amp = DimArray(values(annual_max .- annual_min), (Ti(timestamp(annual_min)), dims(T,Z)))
    zaa = mapslices(T_amp, dims=Z) do A
        i_lo = findfirst(A .< upreferred(threshold))
        if !isnothing(i_lo) && i_lo > 1
            i_hi = i_lo - 1 # one grid cell higher
            z₁ = dims(T, Z)[i_hi]
            z₂ = dims(T, Z)[i_lo]
            A₁ = upreferred(A[i_hi])
            A₂ = upreferred(A[i_lo])
            return solve_depth_linear(threshold, A₁, A₂, z₁, z₂)
        elseif i_lo == 1
            return dims(T,Z)[1]
        else
            return Inf*u"m"
        end
    end
    return rebuild(T, dropdims(zaa, dims=Z), (Ti(timestamp(annual_max)),))
end
"""
    thawdepth(T::AbstractDimArray{<:TempQuantity}; Tmelt=0.0u"°C")

Computes sub-grid thaw depth (a.k.a freezing front) from temperature at all time steps.
The sub-grid depth of the zero degree isotherm is determined by linearly interpolating between
grid cell temperatures in `T`. Note that `T` is assumed to have units °C and dimensions `Ti` (time)
and `Z` (depth) in any order.
"""
function thawdepth(T::AbstractDimArray{<:TempQuantity}; Tmelt=0.0u"°C")
    _check_arr_dims(T)
    T = permutedims(T, (Z,Ti))
    td = mapslices(T, dims=Z) do Tₜ
        i_frozen = findfirst(Tₜ .<= Tmelt)
        i_thawed = findfirst(Tₜ .> Tmelt)
        tdₜ = if !isnothing(i_thawed) && !isnothing(i_frozen) && i_thawed < i_frozen
            z₁ = dims(T, Z)[i_thawed]
            z₂ = dims(T, Z)[i_frozen]
            T₁ = upreferred(Tₜ[i_thawed])
            T₂ = upreferred(Tₜ[i_frozen])
            Tm = upreferred(Tmelt)
            return solve_depth_linear(Tm, T₁, T₂, z₁, z₂)
        elseif isnothing(i_frozen)
            dims(T, Z)[end]
        else
            dims(T, Z)[1]
        end
        return tdₜ
    end
    return rebuild(T, dropdims(td, dims=Z), (dims(T,Ti),))
end
"""
    active_layer_thickness(T::AbstractDimArray{<:TempQuantity})

Computes active layer thickness annually. The active layer thickness is defined here as the maximum thaw
depth throughout the calendar year. Assumes `T` to have dimensions `Ti` (time) and `Z` (depth) in any order.
"""
function active_layer_thickness(T::AbstractDimArray{<:TempQuantity}; Tmelt=0.0u"°C")
    _check_arr_dims(T)
    thaw_depth = thawdepth(T; Tmelt)
    Tarr = TimeArray(dims(T, Ti).val, thaw_depth)
    annual_max = collapse(Tarr, year, first, maximum)
    return rebuild(T, values(annual_max), (Ti(timestamp(annual_max)),))
end
"""
    permafrosttable(T::AbstractDimArray{<:TempQuantity})

Computes depth of permafrost table for all years, i.e. the closest depth to the surface at which
the maximum annual temperature is strictly less than `Tmelt`. Assumes `T` to have dimensions `Ti` (time) and
`Z` (depth) in any order.
"""
function permafrosttable(T::AbstractDimArray{<:TempQuantity}; Tmelt=0.0u"°C")
    _check_arr_dims(T)
    Tarr = TimeArray(dims(T, Ti).val, permutedims(T, (Ti,Z)))
    annual_max = collapse(Tarr, year, first, maximum)
    # find argmax (first one/"true" bit) of annual_max < 0
    z_inds = vec(argmax(values(annual_max) .< Tmelt, dims=2))
    return rebuild(T, collect(dims(T, Z)[[i[2] for i in z_inds]]), (Ti(timestamp(annual_max)),))
end
"""
    permafrostbase(T::AbstractDimArray{<:TempQuantity})

Computes depth of permafrost base for all years, i.e. the closest depth to the "bottom" at which
the maximum annual temperature is strictly `< 0°C`. Assumes `T` to have dimensions `Ti` (time) and
`Z` (depth) in any order.
"""
function permafrostbase(T::AbstractDimArray{<:TempQuantity}; Tmelt=0.0u"°C")
    _check_arr_dims(T)
    Tarr = TimeArray(dims(T, Ti).val, permutedims(T, (Ti,Z)))
    annual_max = collapse(Tarr, year, first, maximum)
    # find argmax (first one/"true" bit) of annual_max < 0
    z_inds = vec(size(Tarr,2) .- [i[2] for i in argmax(reverse(values(annual_max) .< Tmelt, dims=2), dims=2)])
    return rebuild(T, collect(dims(T, Z)[[i for i in z_inds]]), (Ti(timestamp(annual_max)),))
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
