"""
    TemperatureProfileObservable(name::Symbol, zs::AbstractVector, t0::DateTime, tsave::AbstractVector{DateTime}; kwargs...)

Special case of `TileStateObservable` for temperature profiles.
"""
function TemperatureProfileObservable(name::Symbol, zs::AbstractVector, tspan::NTuple{2,DateTime}, p::Period; kwargs...)
    SimulatorObservable(
        name,
        integrator -> collect(map(getvar(Val{:T}(), Tile(integrator), integrator.u; interp=true), ustrip.(zs))),
        tspan[1],
        tspan[1]+p:p:tspan[2],
        (Z(zs),);
        kwargs...
    )
end

function ActiveLayerThicknessObservable(name::Symbol, tspan::NTuple{2,DateTime}; samplerate=Hour(12), kwargs...)
    @assert Date(tspan[2]) - Date(tspan[1]) >= Day(365) "Active layer thickness requires tspan of >= 1 year"
    return SimulatorObservable(
        name,
        sample_thawdepth,
        tspan[1]+samplerate,
        tspan[1]+Year(1):Year(1):tspan[2],
        (1,);
        reducer=maximum,
        samplerate,
        kwargs...
    )
end

function LayerVarObservable(name::Symbol, layername::Symbol, varname::Symbol, grid::Grid, tspan::NTuple{2,DateTime}, p::Period=Day(1); kwargs...)
    SimulatorObservable(
        name,
        integrator -> getproperty(getproperty(getstate(integrator), layername), varname),
        tspan[1],
        tspan[1]+p:p:tspan[2],
        (Z(collect(grid)),);
        kwargs...
    )
end

function sample_thawdepth(integrator)
    state = getstate(integrator)
    tile = Tile(integrator)
    u = integrator.u
    zs = cells(state.grid)
    T_var = getvar(Val{:T}(), tile, u; interp=false).*u"°C"
    T = DimArray(reshape(T_var,:,1), (Z(zs),Ti([convert_t(integrator.t)])))
    td = Diagnostics.thawdepth(T)
    return Vector(td[:,1])
end
