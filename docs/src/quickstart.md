## Quick start

After [installing](installation.md) `CryoGrid.jl`, you can get started right away with a simple soil heat model. The [`Presets`](@ref) module (aliased `CryoGrid.Presets`) provides pre-specified model configurations that can be obtained with a single function call. It is also possible to modify the soil and initial temperature profiles via `SoilLayerConfig`; here `SamoylovDefault` is used.

Single layer heat conduction model with free water freeze curve and air temperature upper boundary condition:

```julia
using CryoGrid
using Plots

# load provided forcing data from Samoylov;
# The forcing file will be automatically downloaded to the input/ folder if not already present.
forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044, :Tair => u"°C");
# use air temperature as upper boundary forcing
tair = TimeSeriesForcing(ustrip.(forcings.data.Tair), forcings.timestamps, :Tair);
# get preset soil and initial temperature profile for Samoylov
soilprofile, tempprofile = CryoGrid.Presets.SamoylovDefault
# basic 1-layer heat conduction model (defaults to free water freezing scheme)
model = CryoGrid.Presets.SoilHeatTile(TemperatureGradient(tair), soilprofile)
# define time span (1 year)
tspan = (DateTime(2010,10,30),DateTime(2011,10,30))
p = parameters(model)
initT = initializer(:T, tempprofile)
u0 = initialcondition!(model, tspan, p, initT)
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(model,u0,tspan,p,savevars=(:T,))
# solve discretized system, saving every 6 hours;
# Trapezoid on a discretized PDE is analogous to the well known Crank-Nicolson method.
out = @time solve(prob, Trapezoid(), saveat=6*3600.0, progress=true) |> CryoGridOutput;
zs = [1.0,5,10,20,30,50,100,500,1000]u"cm"
cg = Plots.cgrad(:copper,rev=true)
plot(out.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false)
```
![Ts_output_freew](res/Ts_H_tair_freeW_2010-2011.png)

Alternatively, we can use a Dall'Amico freeze curve:

```julia
model = CryoGrid.Presets.SoilHeatTile(TemperatureGradient(tair), soilprofile, freezecurve=SFCC(DallAmico()))
# Set-up parameters
p = parameters(model)
tspan = (DateTime(2010,10,30),DateTime(2011,10,30))
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(model,u0,tspan,p,savevars=(:T,))
# stiff solvers don't work well with Dall'Amico due to the ill-conditioned Jacobian;
# We can just forward Euler instead.
out = @time solve(prob, Euler(), dt=120.0, saveat=6*3600.0, progress=true) |> CryoGridOutput;
plot(out.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false)
```

Note that `SoilHeat` uses energy as the state variable by default. To use temperature as the state variable instead:

```julia
# :T is the variable name for temperature, :H represents enthalpy/energy.
# This is used in the specification of the Heat process type.
# While this will work with any freeze curve, here we use Westermann (2011) as an example.
model = CryoGrid.Presets.SoilHeatTile(:T, TemperatureGradient(tair), soilprofile, freezecurve=SFCC(Westermann()))
```