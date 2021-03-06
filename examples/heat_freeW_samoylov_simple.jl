using CryoGrid
using Plots

forcings = loadforcings(CryoGrid.Models.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044, :Tair => u"°C");
# use air temperature as upper boundary forcing;
tair = TimeSeriesForcing(ustrip.(forcings.data.Tair), forcings.timestamps, :Tair);
# basic 1-layer heat conduction model (defaults to free water freezing scheme)
grid = CryoGrid.Models.DefaultGrid_5cm
model = CryoGrid.Models.SoilHeat(:H, TemperatureGradient(tair), CryoGrid.Models.SamoylovDefault; grid=grid)
# define time span
tspan = (DateTime(2010,10,30),DateTime(2011,10,30))
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(model,tspan)
# solve with Crank-Nicolson (Trapezoid) and construct CryoGridOutput from solution
out = @time solve(prob, Trapezoid(), abstol=1e-4, reltol=1e-4, saveat=24*3600.0, progress=true) |> CryoGridOutput;
# Plot it!
zs = [1:10...,20:10:100...]
cg = Plots.cgrad(:copper,rev=true);
plot(out.soil.H[Z(zs)], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Enthalpy", leg=false, dpi=150)
plot(out.soil.T[Z(zs)], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)
