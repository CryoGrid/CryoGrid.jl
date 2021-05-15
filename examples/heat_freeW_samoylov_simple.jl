using CryoGrid
using Plots

forcings = loadforcings(CryoGridModels.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044, :Tair => u"°C");
# use air temperature as upper boundary forcing;
# convert Tair to Kelvin
tair = TimeSeriesForcing(ustrip.(u"K", forcings.data.Tair), forcings.timestamps, :Tair);
# basic 1-layer heat conduction model (defaults to free water freezing scheme)
grid = CryoGridModels.DefaultGrid_5cm
model = CryoGridModels.SoilHeat(:H, TemperatureGradient(tair), CryoGridModels.SamoylovDefault; grid=grid)
# define time span
tspan = (DateTime(2010,10,30),DateTime(2011,10,30))
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(model,tspan)
# solve with Crank-Nicolson (Trapezoid) and construct CryoGridOutput from solution
out = @time solve(prob, Trapezoid(), abstol=1e-4, reltol=1e-4, saveat=6*3600.0, progress=true) |> CryoGridOutput;
# Plot it!
zs = [1:10...,20:10:100...]
cg = Plots.cgrad(:copper,rev=true);
plot(out.soil.H[Z(zs)], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Enthalpy", leg=false, dpi=150)
plot(out.soil.T[Z(zs)] .- 273.15, color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)
