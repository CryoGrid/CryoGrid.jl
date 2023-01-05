using CryoGrid
using Plots

# load provided forcing data from Samoylov;
# The forcing file will be automatically downloaded to the input/ folder if not already present.
forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044, :Tair => u"°C");
# use air temperature as upper boundary forcing
tair = TimeSeriesForcing(forcings.data.Tair, forcings.timestamps, :Tair);
# get preset soil and initial temperature profile for Samoylov
soilprofile, tempprofile = CryoGrid.Presets.SamoylovDefault
initT = initializer(:T, tempprofile)
# basic 1-layer heat conduction model (defaults to free water freezing scheme)
tile = CryoGrid.Presets.SoilHeatTile(TemperatureGradient(tair), soilprofile, initT)
# define time span (1 year)
tspan = (DateTime(2010,10,30),DateTime(2011,10,30))
u0, du0 = initialcondition!(tile, tspan)
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(tile, u0, tspan, savevars=(:T,))
# solve discretized system, saving every 6 hours;
# Trapezoid on a discretized PDE is analogous to the well known Crank-Nicolson method.
out = @time solve(prob, Trapezoid(), saveat=6*3600.0, progress=true) |> CryoGridOutput;
zs = [1.0,5,10,20,30,50,100,500,1000]u"cm"
cg = Plots.cgrad(:copper,rev=true)
plot(out.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false)
savefig("res/Ts_H_tair_freeW_2010-2011.png")

tile = CryoGrid.Presets.SoilHeatTile(TemperatureGradient(tair), soilprofile, freezecurve=SFCC(DallAmico()))
# Set-up parameters
p = parameters(tile)
tspan = (DateTime(2010,10,30),DateTime(2011,10,30))
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(tile,u0,tspan,p,savevars=(:T,))
# stiff solvers don't work well with Dall'Amico due to the ill-conditioned Jacobian;
# We can just forward Euler instead.
out = @time solve(prob, Euler(), dt=120.0, saveat=6*3600.0, progress=true) |> CryoGridOutput;
plot(out.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false)
savefig("res/Ts_H_tair_vg_2010-2011.png")
