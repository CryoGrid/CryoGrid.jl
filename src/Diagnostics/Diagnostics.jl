module Diagnostics

using CryoGrid
using CryoGrid.Numerics
using CryoGrid.Strat
using CryoGrid.Utils

using DimensionalData
using Statistics
using TimeSeries
using Unitful

_check_arr_dims(T::AbstractDimArray) = @assert hasdim(T, Z) && hasdim(T, Ti) "Tay must have depth (Z) and time (Ti) dimensions"
# solve Tm-T₁ = m*(z - z₁) for z where m = (T₂ - T₁)/(z₂ - z₁)
solve_depth_linear(T, T₁, T₂, z₁, z₂) = (T-T₁)*(z₂ - z₁)/(T₂ - T₁) + z₁

export integrate, computejac
include("numerics.jl")

export zero_annual_amplitude, permafrosttable, permafrostbase, thawdepth, active_layer_thickness, mean_annual_ground_temperature
include("permafrost.jl")

export spinup
include("spinup.jl")

end