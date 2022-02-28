"""
Pre-built CryoGrid configurations for rapid prototyping.
"""
module Presets

using CryoGrid
using CryoGrid.InputOutput: Resource
using CryoGrid.Numerics

using Statistics

export SoilHeatColumn, SamoylovDefault

include("presetgrids.jl")

"""
    SoilHeatColumn([heatvar=:H], upperbc::BoundaryProcess, soilprofile::Profile, init::Numerics.VarInitializer; grid::Grid=DefaultGrid, freezecurve::F=FreeWater()) where {F<:FreezeCurve}

Builds a simple one-layer soil/heat-conduction model with the given grid and configuration. Uses the "free water" freeze curve by default,
but this can be changed via the `freezecurve` parameter. For example, to use the Dall'Amico freeze curve, set `freezecurve=SFCC(DallAmico())`.
"""
function SoilHeatColumn(heatvar, upperbc::BoundaryProcess, soilprofile::Profile, init::Numerics.VarInitializer;
    grid::Grid=DefaultGrid_5cm, freezecurve::F=FreeWater(), chunksize=nothing) where {F<:FreezeCurve}
    strat = Stratigraphy(
        -2.0u"m" => top(upperbc),
        Tuple(z => subsurface(Symbol(:soil,i), Soil(para=para), Heat(heatvar, freezecurve=freezecurve)) for (i,(z,para)) in enumerate(soilprofile)),
        1000.0u"m" => bottom(GeothermalHeatFlux(0.053u"J/s/m^2"))
    )
    Tile(strat, grid, init, chunksize=chunksize)
end
SoilHeatColumn(upperbc::BoundaryProcess, soilprofile::Profile, init::Numerics.VarInitializer; grid::Grid=DefaultGrid_2cm, freezecurve::F=FreeWater()) where {F<:FreezeCurve} = SoilHeatColumn(:H, upperbc, soilprofile, init; grid=grid, freezecurve=freezecurve)

Forcings = (
    Samoylov_ERA5_fitted_daily_1979_2020 = Resource("samoylov_era5_fitted_daily_1979-2020", "json", "https://nextcloud.awi.de/s/DSx3BWsfjCdy9MF"),
    Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044 = Resource("Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044", "json", "https://nextcloud.awi.de/s/dp74KC2ceQKaG43"),
    Samoylov_ERA_MkL3_CCSM4_long_term = Resource("Samoylov_ERA_MkL3_CCSM4_long_term", "json", "https://nextcloud.awi.de/s/gyoMTy9jpk2pMxL"),
)
Parameters = (
    # Faroux et al. doi:10.1109/IGARSS.2007.4422971
    EcoCLimMap_ULC_126_72 = Resource("EcoCLimMap_ULC_126_72", "json", "https://nextcloud.awi.de/s/nWiJr5pBoqFtw7p/download")
)

const SamoylovDefault = (
    soilprofile = SoilProfile(
        0.0u"m" => CharacteristicFractions(xic=0.0,por=0.80,sat=1.0,org=0.75), #(θw=0.80,θm=0.05,θo=0.15,ϕ=0.80),
        0.1u"m" => CharacteristicFractions(xic=0.0,por=0.80,sat=1.0,org=0.25), #(θw=0.80,θm=0.15,θo=0.05,ϕ=0.80),
        0.4u"m" => CharacteristicFractions(xic=0.30,por=0.55,sat=1.0,org=0.25), #(θw=0.80,θm=0.15,θo=0.05,ϕ=0.55),
        3.0u"m" => CharacteristicFractions(xic=0.0,por=0.50,sat=1.0,org=0.0), #(θw=0.50,θm=0.50,θo=0.0,ϕ=0.50),
        10.0u"m" => CharacteristicFractions(xic=0.0,por=0.30,sat=1.0,org=0.0), #(θw=0.30,θm=0.70,θo=0.0,ϕ=0.30),
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