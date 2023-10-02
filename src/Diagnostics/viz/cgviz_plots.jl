function plot_at_depths(var::Symbol, out::CryoGridOutput, depths::AbstractVector; cmap=Plots.cgrad(:copper, rev=true), kwargs...)
    data = getproperty(out, var)
    Plots.plot(data[Z(Near(depths))]; color=cmap[LinRange(0.0,1.0,length(depths))]', leg=false, kwargs...)
end
