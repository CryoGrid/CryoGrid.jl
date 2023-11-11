"""
    WaterTableInitializer <: VarInitializer{:sat}

Simple, piecewise constant initializer for saturation state that takes a surface-level
saturation `sat0` and water table depth `z_tab` and produces a two-segment, piecewise
constant profile with the saturation level set to `(sat0 + 1.0) / 2` from the halfway
point down to the water table.
"""
Base.@kwdef struct WaterTableInitializer{Tsat,Tz0,Tzp,TI} <: VarInitializer{:sat}
    sat0::Tsat = 0.5
    z0::Tz0 = 0.0u"m"
    z_tab::Tzp = 10.0u"m"
    interp::TI = Interp.Constant()
end

(init::WaterTableInitializer)(::SubSurface, state) = init(state.sat, state.grid)

function (init::WaterTableInitializer)(sat::AbstractVector, grid::Grid)
    z0 = ustrip(init.z0)
    z_tab = ustrip(init.z_tab)
    sat0 = ustrip(init.sat0)
    f = Interp.extrapolate(
        Interp.interpolate(
            ([z0, (z_tab - z0)/2, z_tab],),
            [sat0, (sat0 + 1.0)/2, 1.0],
            Interp.Gridded(init.interp)
        ),
        Interp.Flat(),
    )
    sat .= f.(cells(grid))
end
