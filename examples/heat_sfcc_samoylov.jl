using CryoGrid
using Plots

forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044, :Tair => u"°C");
# use air temperature as upper boundary forcing;
tair = TimeSeriesForcing(forcings.data.Tair, forcings.timestamps, :Tair);
# "simple" heat conduction model w/ 5 cm grid spacing
grid = CryoGrid.Presets.DefaultGrid_5cm
soilprofile, tempprofile = CryoGrid.Presets.SamoylovDefault
initT = initializer(:T, tempprofile)
sfcc = PainterKarra(swrc=VanGenuchten(α=0.1, n=2.0))
tile = CryoGrid.Presets.SoilHeatTile(:H, TemperatureGradient(tair), GeothermalHeatFlux(0.053u"W/m^2"), soilprofile, initT; grid=grid, freezecurve=sfcc)
# define time span
tspan = (DateTime(2010,11,30),DateTime(2011,11,30))
u0, du0 = initialcondition!(tile, tspan)
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(tile, u0, tspan, savevars=(:T,))
@info "Running model"
# solve with adaptive strong stability preserving Runge-Kutta method and construct CryoGridOutput from solution
out = @time solve(prob, Euler(), dt=600.0, saveat=3*3600.0, progress=true) |> CryoGridOutput;
# Plot it!
zs = [5,10,15,20,25,30,40,50,100,500,1000]u"cm"
cg = Plots.cgrad(:copper,rev=true);
plot(out.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)
