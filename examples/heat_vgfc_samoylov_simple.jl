using CryoGrid
using Plots

forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044, :Tair => u"°C");
# use air temperature as upper boundary forcing;
tair = TimeSeriesForcing(ustrip.(forcings.data.Tair), forcings.timestamps, :Tair);
# "simple" heat conduction model w/ 5 cm grid spacing
grid = CryoGrid.Presets.DefaultGrid_5cm
soilprofile, tempprofile = CryoGrid.Presets.SamoylovDefault
initT = initializer(:T, tempprofile)
model = CryoGrid.Presets.SoilHeatColumn(:H, TemperatureGradient(tair), soilprofile, initT; grid=grid, freezecurve=SFCC(DallAmico()))
# define time span
tspan = (DateTime(2010,10,30),DateTime(2011,10,30))
p = parameters(model)
u0, du0 = initialcondition!(model, tspan, p)
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(model,u0,tspan,p,savevars=(:T,))
# solve with forward Euler (fixed 10 minute timestep) and construct CryoGridOutput from solution
out = @time solve(prob, Euler(), dt=10*60.0, saveat=24*3600.0, progress=true) |> CryoGridOutput;
# Plot it!
zs = [1:10...,20:10:100...]
cg = Plots.cgrad(:copper,rev=true);
plot(out.H[Z(zs)], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Enthalpy", leg=false, dpi=150)
plot(out.T[Z(zs)], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)