"""
Pre-built CryoGrid configurations for rapid prototyping.
"""
module Models

using ..CryoGrid

include("presetgrids.jl")

Forcings = (
    Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044 = CryoGrid.Resource("Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044", "json", "https://nextcloud.awi.de/s/F98s5WEo9xMPod7/download"),
    Samoylov_ERA_MkL3_CCSM4_long_term = CryoGrid.Resource("Samoylov_ERA_MkL3_CCSM4_long_term", "json", "https://nextcloud.awi.de/s/RSqBtp5sPwkCf45/download"),
)

Parameters = (
    # Faroux et al. doi:10.1109/IGARSS.2007.4422971
    EcoCLimMap_ULC_126_72 = CryoGrid.Resource("EcoCLimMap_ULC_126_72", "json", "https://nextcloud.awi.de/s/nWiJr5pBoqFtw7p/download")
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
    SoilProfile(
        0.0u"m" => (0.80,0.0,0.05,0.15,0.80),
        0.1u"m" => (0.80,0.0,0.15,0.05,0.80),
        0.4u"m" => (0.80,0.0,0.15,0.05,0.55),
        3.0u"m" => (0.50,0.0,0.50,0.0,0.50),
        10.0u"m" => (0.30,0.0,0.70,0.0,0.30),
    ),
    TempProfile(
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
    SoilHeat([heatvar=:H], upperbc::BoundaryProcess{Heat}, soilconfig::SoilLayerConfig; grid::Grid=DefaultGrid, freezecurve::F=FreeWater()) where {F<:FreezeCurve}

Builds a simple one-layer soil/heat-conduction model with the given grid and configuration. Uses the "free water" freeze curve by default,
but this can be changed via the `freezecurve` parameter. For example, to use the van Genuchten freeze curve, set `freezecurve=SFCC(VanGenuchten())`.
"""
function SoilHeat(heatvar, upperbc::BoundaryProcess{Heat}, soilconfig::SoilLayerConfig;
    grid::Grid=DefaultGrid_5cm, freezecurve::F=FreeWater(), chunk_size=nothing) where {F<:FreezeCurve}
    strat = Stratigraphy(
        -2.0u"m" => Top(upperbc),
        0.0u"m" => Ground(:soil, Soil{Sand}(soilconfig.soilprofile), Heat{heatvar}(soilconfig.tempprofile, freezecurve=freezecurve)),
        1000.0u"m" => Bottom(GeothermalHeatFlux(0.053u"J/s/m^2"))
    )
    model = CryoGridSetup(strat,grid,chunk_size=chunk_size)
end
SoilHeat(upperbc::BoundaryProcess{Heat}, soilconfig::SoilLayerConfig; grid::Grid=DefaultGrid_2cm, freezecurve::F=FreeWater()) where {F<:FreezeCurve} =
    SoilHeat(:H, upperbc, soilconfig; grid=grid, freezecurve=freezecurve)
end