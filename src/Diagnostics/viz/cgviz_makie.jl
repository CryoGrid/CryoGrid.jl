function plot_at_depths(var::Symbol, out::CryoGridOutput, depths; kwargs...)
    fig = Makie.Figure()
    plot_at_depths!(fig[1,1], var, out, depths; kwargs...)
    return fig
end
function plot_at_depths!(fig, var::Symbol, out::CryoGridOutput, depths; axis_kwargs=(;), ylabel="", title="", xticks=nothing, cmap::Symbol=:copper, kwargs...)
    data = getproperty(out, var)
    subdata = data[Z(Near(depths))]
    ts = collect(dims(data, Ti))
    xticks = if isnothing(xticks)
        # auto choose x tick interval as N / log(N)
        intrv = Int(floor(length(ts) / ceil(log(length(ts)))))
        (1:intrv:length(ts), Dates.format.(ts[1:intrv:length(ts)], "YYYY-mm-dd"))
    else
        xticks
    end
    ax = Makie.Axis(fig; title, ylabel, xticks, xticklabelrotation = π/4, axis_kwargs...)
    cm = Makie.cgrad(cmap, size(subdata, Z), rev=true, categorical=true)
    Makie.series!(1:length(ts), ustrip.(subdata.data), color=collect(cm), kwargs...)
    return ax
end
