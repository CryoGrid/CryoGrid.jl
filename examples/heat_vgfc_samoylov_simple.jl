using CryoGrid
using Plots

forcings = loadforcings(CryoGrid.Models.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044, :Tair => u"°C");
# use air temperature as upper boundary forcing;
tair = TimeSeriesForcing(ustrip.(forcings.data.Tair), forcings.timestamps, :Tair);
# basic 1-layer heat conduction model (defaults to free water freezing scheme)
grid = CryoGrid.Models.DefaultGrid_5cm
model = CryoGrid.Models.SoilHeat(:H, TemperatureGradient(tair), CryoGrid.Models.SamoylovDefault;
    grid=grid, freezecurve=SFCC(DallAmico()))
# define time span
tspan = (DateTime(2010,10,30),DateTime(2011,10,30))
p = copy(model.pproto)
p.soil.α .= 4.0
p.soil.n .= 2.0
p.soil.Tₘ .= 0.0
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(model,tspan,p)
# solve with forward Euler (fixed 10 minute timestep) and construct CryoGridOutput from solution
out = @time solve(prob, Euler(), dt=10*60.0, saveat=24*3600.0, progress=true) |> CryoGridOutput;
# Plot it!
zs = [1:10...,20:10:100...]
cg = Plots.cgrad(:copper,rev=true);
plot(out.soil.H[Z(zs)], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Enthalpy", leg=false, dpi=150)
plot(out.soil.T[Z(zs)], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)