"""
Pre-built CryoGrid configurations for rapid prototyping.
"""
module Presets

using CryoGrid
using CryoGrid.InputOutput: Resource
using CryoGrid.Numerics
using CryoGrid.Utils

# physics modules
using CryoGrid.Heat
using CryoGrid.Hydrology
using CryoGrid.Soils

using Statistics

export SoilHeatTile, SamoylovDefault

include("presetgrids.jl")

"""
    SoilHeatTile([heatop=:H], upperbc::BoundaryProcess, soilprofile::Profile, init::CryoGrid.Initializer...; grid::Grid=DefaultGrid_10cm, tile_kwargs...) where {F<:FreezeCurve}

Builds a simple one-layer soil/heat-conduction model with the given grid and configuration.
"""
function SoilHeatTile(heatop, upperbc::BoundaryProcess, lowerbc::BoundaryProcess, soilprofile::Profile, inits::CryoGrid.Initializer...; grid::Grid=DefaultGrid_5cm, tile_kwargs...)
    strat = Stratigraphy(
        grid[1] => Top(upperbc),
        Tuple(d => Ground(para, heat=HeatBalance(heatop), water=nothing) for (i,(d,para)) in enumerate(soilprofile)),
        grid[end] => Bottom(lowerbc)
    )
    return Tile(strat, PresetGrid(grid), inits...; tile_kwargs...)
end
SoilHeatTile(upperbc::BoundaryProcess, lowerbc::BoundaryProcess, soilprofile::Profile, inits::CryoGrid.Initializer...; kwargs...) = SoilHeatTile(:H, upperbc, lowerbc, soilprofile, inits...; kwargs...)

Forcings = (
    Samoylov_ERA5_fitted_daily_1979_2020 = Resource("samoylov_era5_fitted_daily_1979-2020", ForcingFormatJSON{2}(), "https://nextcloud.awi.de/s/ScYAoHzeMzAfpjf/download/samoylov_era5_fitted_daily_1979-2020.json"),
    Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044 = Resource("Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044", ForcingFormatJSON{1}(), "https://nextcloud.awi.de/s/cbeycGQoQpXi3Ei/download/samoylov_ERA_obs_fitted_1979_2014_spinup_extended2044.json"),
    Samoylov_ERA_MkL3_CCSM4_long_term = Resource("Samoylov_ERA_MkL3_CCSM4_long_term", ForcingFormatJSON{1}(), "https://nextcloud.awi.de/s/45ax9AsTACxL25Q/download/FORCING_ULC_126_72.json"),
    Bayelva_ERA5_fitted_daily_1979_2020 = Resource("bayelva_era5_fitted_daily_1979-2020", ForcingFormatJSON{2}(), "https://nextcloud.awi.de/s/5AdbRMYKneCHgx4/download/bayelva_era5_fitted_daily_1979-2020.json")
)
Parameters = (
    # Faroux et al. doi:10.1109/IGARSS.2007.4422971
    EcoCLimMap_ULC_126_72 = Resource("EcoCLimMap_ULC_126_72", ParamsJSON{1}(), "https://nextcloud.awi.de/s/7F65JET9TzdosMD/download/PARA_ULC_126_72.json")
)

const SamoylovDefault = (
    soilprofile = SoilProfile(
        0.0u"m" => SimpleSoil(por=0.80,sat=1.0,org=0.75,freezecurve=PainterKarra(swrc=VanGenuchten("organic"))),
        0.1u"m" => SimpleSoil(por=0.80,sat=1.0,org=0.25,freezecurve=PainterKarra(swrc=VanGenuchten("organic"))),
        0.4u"m" => SimpleSoil(por=0.55,sat=1.0,org=0.25,freezecurve=PainterKarra(swrc=VanGenuchten("sandy loam"))),
        3.0u"m" => SimpleSoil(por=0.50,sat=1.0,org=0.0,freezecurve=PainterKarra(swrc=VanGenuchten("sandy loam"))),
        10.0u"m" => SimpleSoil(por=0.30,sat=1.0,org=0.0,freezecurve=PainterKarra(swrc=VanGenuchten("sandy loam"))),
        30.0u"m" => SimpleSoil(por=0.05,sat=1.0,org=0.0,freezecurve=PainterKarra(swrc=VanGenuchten("sand"))),
    ),
    tempprofile = TemperatureProfile(
        0.0u"m" => -1.0u"°C",
        0.2u"m" => -1.0u"°C",
        2.0u"m" => -1.0u"°C",
        5.0u"m" => -3.0u"°C",
        9.5u"m" => -6.0u"°C",
        25.0u"m" => -9.0u"°C",
        100.0u"m" => -9.0u"°C",
        1000.0u"m" => 10.2u"°C",
    )
)

end