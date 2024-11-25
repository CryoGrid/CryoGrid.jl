"""
    SoilHeatTile([heatop=:H], upperbc::BoundaryProcess, soilprofile::Profile, init::CryoGrid.Initializer...; grid::Grid=DefaultGrid_10cm, tile_kwargs...) where {F<:FreezeCurve}

Builds a simple one-layer soil/heat-conduction model with the given grid and configuration.
"""
function SoilHeatTile(
    heatop,
    upperbc::BoundaryProcess,
    lowerbc::BoundaryProcess,
    soilprofile::Profile,
    inputs::InputProvider,
    inits::CryoGrid.Initializer...;
    grid::Grid=DefaultGrid_5cm,
    tile_kwargs...
)
    strat = Stratigraphy(
        grid[1] => Top(upperbc),
        Tuple(d => Ground(para, heat=HeatBalance(heatop), water=nothing) for (i,(d,para)) in enumerate(soilprofile)),
        grid[end] => Bottom(lowerbc)
    )
    return Tile(strat, PresetGrid(grid), inputs, inits...; tile_kwargs...)
end
SoilHeatTile(upperbc::BoundaryProcess, lowerbc::BoundaryProcess, soilprofile::Profile, inputs::InputProvider, inits::CryoGrid.Initializer...; kwargs...) = SoilHeatTile(:H, upperbc, lowerbc, soilprofile, inputs, inits...; kwargs...)

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
