"""
Pre-built CryoGrid configurations for rapid prototyping.
"""
module Presets

using CryoGrid
using CryoGrid.InputOutput: Resource
using CryoGrid.Numerics

using Statistics

export SoilHeatTile, SamoylovDefault

include("presetgrids.jl")

"""
    SoilHeatTile([heatop=:H], upperbc::BoundaryProcess, soilprofile::Profile, init::VarInitializer; grid::Grid=DefaultGrid_10cm, freezecurve::F=FreeWater()) where {F<:FreezeCurve}

Builds a simple one-layer soil/heat-conduction model with the given grid and configuration. Uses the "free water" freeze curve by default,
but this can be changed via the `freezecurve` parameter. For example, to use the Dall'Amico freeze curve, set `freezecurve=SFCC(DallAmico())`.
"""
function SoilHeatTile(heatop, upperbc::BoundaryProcess, lowerbc::BoundaryProcess, soilprofile::Profile, init::VarInitializer; grid::Grid=DefaultGrid_5cm, freezecurve::F=FreeWater(), chunk_size=nothing) where {F<:FreezeCurve}
    strat = Stratigraphy(
        grid[1] => Top(upperbc),
        Tuple(knot.depth => Symbol(:soil,i) => HomogeneousSoil(knot.value, heat=HeatBalance(heatop, freezecurve=freezecurve)) for (i,knot) in enumerate(soilprofile)),
        grid[end] => Bottom(lowerbc)
    )
    return Tile(strat, PresetGrid(grid), init, chunk_size=chunk_size)
end
SoilHeatTile(upperbc::BoundaryProcess, lowerbc::BoundaryProcess, soilprofile::Profile, init::VarInitializer; grid::Grid=DefaultGrid_5cm, freezecurve::F=FreeWater()) where {F<:FreezeCurve} = SoilHeatTile(:H, upperbc, lowerbc, soilprofile, init; grid=grid, freezecurve=freezecurve)

Forcings = (
    Samoylov_ERA5_fitted_daily_1979_2020 = Resource("samoylov_era5_fitted_daily_1979-2020", ForcingFormatJSON{2}(), "https://nextcloud.awi.de/s/HABRPPxekqgHsSC/download/samoylov_era5_fitted_daily_1979-2020.json"),
    Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044 = Resource("Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044", ForcingFormatJSON{1}(), "https://nextcloud.awi.de/s/kcjZFstkPcAdQmM/download/samoylov_ERA_obs_fitted_1979_2014_spinup_extended2044.json"),
    Samoylov_ERA_MkL3_CCSM4_long_term = Resource("Samoylov_ERA_MkL3_CCSM4_long_term", ForcingFormatJSON{1}(), "https://nextcloud.awi.de/s/ZiNr3L5XteTe2bN/download/FORCING_ULC_126_72.json"),
    Bayelva_ERA5_fitted_daily_1979_2020 = Resource("bayelva_era5_fitted_daily_1979-2020", ForcingFormatJSON{2}(), "https://nextcloud.awi.de/s/F2i5DBfKCCPKE7m/download/bayelva_era5_fitted_daily_1979-2020.json")
)
Parameters = (
    # Faroux et al. doi:10.1109/IGARSS.2007.4422971
    EcoCLimMap_ULC_126_72 = Resource("EcoCLimMap_ULC_126_72", ParamsJSON{1}(), "https://nextcloud.awi.de/s/7F65JET9TzdosMD/download/PARA_ULC_126_72.json")
)

const SamoylovDefault = (
    soilprofile = SoilProfile(
        0.0u"m" => MineralOrganic(por=0.80,sat=1.0,org=0.75), #(θwi=0.80,θm=0.05,θo=0.15,ϕ=0.80),
        0.1u"m" => MineralOrganic(por=0.80,sat=1.0,org=0.25), #(θwi=0.80,θm=0.15,θo=0.05,ϕ=0.80),
        0.4u"m" => MineralOrganic(por=0.55,sat=1.0,org=0.25), #(θwi=0.80,θm=0.15,θo=0.05,ϕ=0.55),
        3.0u"m" => MineralOrganic(por=0.50,sat=1.0,org=0.0), #(θwi=0.50,θm=0.50,θo=0.0,ϕ=0.50),
        10.0u"m" => MineralOrganic(por=0.30,sat=1.0,org=0.0), #(θwi=0.30,θm=0.70,θo=0.0,ϕ=0.30),
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