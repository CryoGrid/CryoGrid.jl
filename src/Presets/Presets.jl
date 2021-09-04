"""
Pre-built CryoGrid configurations for rapid prototyping.
"""
module Presets

using CryoGrid
using CryoGrid.InputOutput: Resource
using Statistics

include("presetgrids.jl")

Forcings = (
      Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044 = Resource("Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044", "json", "https://nextcloud.awi.de/s/F98s5WEo9xMPod7/download"),
      Samoylov_ERA_MkL3_CCSM4_long_term = Resource("Samoylov_ERA_MkL3_CCSM4_long_term", "json", "https://nextcloud.awi.de/s/RSqBtp5sPwkCf45/download"),
)
Parameters = (
    # Faroux et al. doi:10.1109/IGARSS.2007.4422971
    EcoCLimMap_ULC_126_72 = Resource("EcoCLimMap_ULC_126_72", "json", "https://nextcloud.awi.de/s/nWiJr5pBoqFtw7p/download")
)

"""
    SoilLayerConfig

Helper type for representing site-specific soil layer configuration (e.g. soil and temperature profile).
"""
struct SoilLayerConfig{TSoilProfile,TTempProfile}
    soilprofile::TSoilProfile
    tempprofile::TTempProfile
end

const SamoylovDefault = SoilLayerConfig(
    # soil profile: depth => (total water, liquid water, mineral organic, porosity)
    (
        0.0u"m" => SoilComposition(xic=0.0,por=0.80,sat=1.0,org=0.75), #(θw=0.80,θm=0.05,θo=0.15,ϕ=0.80),
        0.1u"m" => SoilComposition(xic=0.0,por=0.80,sat=1.0,org=0.25), #(θw=0.80,θm=0.15,θo=0.05,ϕ=0.80),
        0.4u"m" => SoilComposition(xic=0.30,por=0.55,sat=1.0,org=0.25), #(θw=0.80,θm=0.15,θo=0.05,ϕ=0.55),
        3.0u"m" => SoilComposition(xic=0.0,por=0.50,sat=1.0,org=0.0), #(θw=0.50,θm=0.50,θo=0.0,ϕ=0.50),
        10.0u"m" => SoilComposition(xic=0.0,por=0.30,sat=1.0,org=0.0), #(θw=0.30,θm=0.70,θo=0.0,ϕ=0.30),
    ),
    TemperatureProfile(
        0.0u"m" => -1.0u"°C",
        2.0u"m" => -1.0u"°C",
        5.0u"m" => -3.0u"°C",
        10.0u"m" => -6.0u"°C",
        25.0u"m" => -9.0u"°C",
        100.0u"m" => -9.0u"°C",
        1000.0u"m" => 10.2u"°C",
    )
)

export SamoylovDefault

"""
    SoilHeat([heatvar=:H], upperbc::BoundaryProcess, soilconfig::SoilLayerConfig; grid::Grid=DefaultGrid, freezecurve::F=FreeWater()) where {F<:FreezeCurve}

Builds a simple one-layer soil/heat-conduction model with the given grid and configuration. Uses the "free water" freeze curve by default,
but this can be changed via the `freezecurve` parameter. For example, to use the Dall'Amico freeze curve, set `freezecurve=SFCC(DallAmico())`.
"""
function SoilHeat(heatvar, upperbc::BoundaryProcess, soilconfig::SoilLayerConfig;
    grid::Grid=DefaultGrid_5cm, freezecurve::F=FreeWater(), chunk_size=nothing) where {F<:FreezeCurve}
    strat = Stratigraphy(
        -2.0u"m" => top(upperbc),
        Tuple(z => subsurface(Symbol(:soil,i), Soil(comp=comp), Heat(sp=heatvar,initialT=soilconfig.tempprofile, freezecurve=freezecurve)) for (i,(z,comp)) in enumerate(soilconfig.soilprofile)),
        1000.0u"m" => bottom(GeothermalHeatFlux(0.053u"J/s/m^2"))
    )
    CryoGridSetup(strat,grid,chunk_size=chunk_size)
end
SoilHeat(upperbc::BoundaryProcess, soilconfig::SoilLayerConfig; grid::Grid=DefaultGrid_2cm, freezecurve::F=FreeWater()) where {F<:FreezeCurve} =
    SoilHeat(:H, upperbc, soilconfig; grid=grid, freezecurve=freezecurve)

end