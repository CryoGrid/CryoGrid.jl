using CryoGrid
using Plots

forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044, :Tair => u"°C");
# use air temperature as upper boundary forcing;
tair = TimeSeriesForcing(ustrip.(forcings.data.Tair), forcings.timestamps, :Tair);
# use default profiles for samoylov
soilprofile, tempprofile = CryoGrid.Presets.SamoylovDefault
# "simple" heat conduction model w/ 5 cm grid spacing (defaults to free water freezing scheme)
grid = CryoGrid.Presets.DefaultGrid_5cm
initT = initializer(:T, tempprofile)
model = CryoGrid.Presets.SoilHeatColumn(:H, TemperatureGradient(tair), soilprofile, initT; grid=grid)
# define time span
tspan = (DateTime(2010,10,30),DateTime(2011,10,30))
p = parameters(model)
u0, du0 = initialcondition!(model, tspan, p)
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(model,u0,tspan,p,savevars=(:T,))
# solve with Crank-Nicolson (Trapezoid) and construct CryoGridOutput from solution
out = @time solve(prob, Trapezoid(), abstol=1e-4, reltol=1e-4, saveat=24*3600.0, progress=true) |> CryoGridOutput;
# Plot it!
zs = [1:10...,20:10:100...]
cg = Plots.cgrad(:copper,rev=true);
plot(out.H[Z(zs)], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Enthalpy", leg=false, dpi=150)
plot(out.T[Z(zs)], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)
