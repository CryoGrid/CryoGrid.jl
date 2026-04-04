"""
    TemperatureProfileObservable(name::Symbol, zs::AbstractVector, t0::DateTime, tsave::AbstractVector{DateTime}; kwargs...)

Special case of `TileStateObservable` for temperature profiles.
"""
function TemperatureProfileObservable(name::Symbol, zs::AbstractVector, tspan::NTuple{2,DateTime}, p::Period; kwargs...)
    SimulatorObservable(
        integrator -> collect(map(getvar(Val{:T}(), Tile(integrator), integrator.u; interp=true), ustrip.(zs))),
        (length(zs),),
        name = name,
        output = TimeSampled(tspan[1], tspan[1]+p:p:tspan[2]),
        time_converter = CryoGrid.convert_t,
        kwargs...
    )
end

function ActiveLayerThicknessObservable(name::Symbol, tspan::NTuple{2,DateTime}; samplerate=Hour(12), kwargs...)
    @assert Date(tspan[2]) - Date(tspan[1]) >= Day(365) "Active layer thickness requires tspan of >= 1 year"
    return SimulatorObservable(
        sample_thawdepth,
        (1,),
        name = name,
        output = TimeSampled(tspan[1]+samplerate, tspan[1]+Year(1):Year(1):tspan[2], samplerate=samplerate),
        time_converter = CryoGrid.convert_t,
        reducer = maximum,
        kwargs...
    )
end

function LayerVarObservable(name::Symbol, layername::Symbol, varname::Symbol, grid::Grid, tspan::NTuple{2,DateTime}, p::Period=Day(1); kwargs...)
    SimulatorObservable(
        integrator -> getproperty(getproperty(getstate(integrator), layername), varname),
        (length(collect(grid)),),
        name = name,
        output = TimeSampled(tspan[1], tspan[1]+p:p:tspan[2]),
        time_converter = CryoGrid.convert_t,
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
