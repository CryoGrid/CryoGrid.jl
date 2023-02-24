using CryoGrid
using Plots

forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044, :Tair => u"°C");
# use air temperature as upper boundary forcing;
tair = TimeSeriesForcing(forcings.data.Tair, forcings.timestamps, :Tair);
# use default profiles for samoylov
soilprofile, tempprofile = CryoGrid.Presets.SamoylovDefault
# "simple" heat conduction model w/ 5 cm grid spacing (defaults to free water freezing scheme)
grid = CryoGrid.Presets.DefaultGrid_5cm
initT = initializer(:T, tempprofile)
tile = CryoGrid.Presets.SoilHeatTile(:H, TemperatureGradient(tair), GeothermalHeatFlux(0.053u"W/m^2"), soilprofile, initT; grid=grid)
# define time span
tspan = (DateTime(2010,10,30),DateTime(2011,10,30))
u0, du0 = initialcondition!(tile, tspan)
# construct CryoGridProblem with tile, initial condition, and timespan;
# we disable the default timestep limiter since we will use an adaptive solver.
prob = CryoGridProblem(tile, u0, tspan, savevars=(:T,), step_limiter=nothing)
@info "Running model"
# solve with Crank-Nicolson (trapezoid method) and construct CryoGridOutput from solution
out = @time solve(prob, Trapezoid(), saveat=24*3600.0, progress=true) |> CryoGridOutput;
# Plot it!
zs = [1,10,20,30,50,100,200,500,1000]u"cm"
cg = Plots.cgrad(:copper,rev=true);
plot(out.H[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Enthalpy", leg=false, dpi=150)
plot(out.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)
