"""
Pre-built CryoGrid configurations for rapid prototyping.
"""
module Models

using ..CryoGrid

const DefaultGrid = Grid(vcat([
            0:0.02:2...,
            2.05:0.05:4.0...,
            4.1:0.1:10...,
            10.2:0.2:20...,
            21:1:30...,
            35:5:50...,
            60:10:100...,
            200:100:1000...
        ]...
    )
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
    GroundHeatOnly(upperbc::BoundaryProcess{Heat}, soilconfig::SoilLayerConfig; grid::Grid=DefaultGrid, freezecurve::F=FreeWater()) where {F<:FreezeCurve}

Builds a simple one-layer soil/heat-conduction model with the given grid and configuration. Uses the "free water" freeze curve by default,
but this can be changed via the `freezecurve` parameter. For example, to use the van Genuchten freeze curve, set `freezecurve=SFCC(VanGenuchten())`.
"""
function GroundHeatOnly(upperbc::BoundaryProcess{Heat}, soilconfig::SoilLayerConfig;
    grid::Grid=DefaultGrid, freezecurve::F=FreeWater(), hcunit=u"J") where {F<:FreezeCurve}
    strat = Stratigraphy(
        -2.0u"m" => Top(upperbc),
        0.0u"m" => Ground(:soil, Soil{Sand}(soilconfig.soilprofile), Heat{hcunit}(soilconfig.tempprofile, freezecurve=freezecurve)),
        1000.0u"m" => Bottom(GeothermalHeatFlux(0.05u"J/s/m^2"))
    )
    model = CryoGridSetup(strat,grid)
end

end