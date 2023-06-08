# method stubs for plotting routines

"""
    plot_at_depths(var::Symbol, out::CryoGridOutput, depths::AbstractArray; kwargs...)

Plots a time series of `var` from `out` at the given `depths`.
"""
function plot_at_depths end

"""
    plot_alt(out::CryoGridOutput; kwargs...)

Plots the active layer thickness from the given `CryoGridOutput`.
"""
function plot_alt end

"""
    plot_temperature_heatmap(out::CryoGridOutput; kwargs...)

Plots a heatmap of temperatures along with the zero degree isotherm.
"""
function plot_temperature_heatmap end
