module Diagnostics

using CryoGrid
using CryoGrid.Numerics
using CryoGrid.Utils

using Dates
using DimensionalData
using Statistics
using TimeSeries
using Unitful
using Requires

_check_arr_dims(T::AbstractDimArray) = @assert hasdim(T, Z) && hasdim(T, Ti) "Tay must have depth (Z) and time (Ti) dimensions"
# solve Tm-T₁ = m*(z - z₁) for z where m = (T₂ - T₁)/(z₂ - z₁)
solve_depth_linear(T, T₁, T₂, z₁, z₂) = (T-T₁)*(z₂ - z₁)/(T₂ - T₁) + z₁

export integrate, computejac
include("numerics.jl")

export zero_annual_amplitude, permafrosttable, permafrostbase, thawdepth, active_layer_thickness, mean_annual_ground_temperature
include("permafrost.jl")

export spinup
include("spinup.jl")

include("viz/cgviz.jl")

function __init__()
    @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" begin
        import .Makie
        include("viz/cgviz_makie.jl")
    end
    @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
        import .Plots
        include("viz/cgviz_plots.jl")
    end
end

end