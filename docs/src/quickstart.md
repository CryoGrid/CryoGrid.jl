## Quick start

After [installing](installation.md) `CryoGrid.jl`, you can get started right away with a simple soil heat model. The [`Presets`](@ref) module (aliased `CryoGrid.Presets`) provides pre-specified model configurations that can be obtained with a single function call. It is also possible to provide a custom soil and temperature profile to `SoilHeatTile`; here the values given in `Presets.SamoylovDefault` is used.

```julia
using CryoGrid
using Plots

# load provided forcing data from Samoylov;
# The forcing file will be automatically downloaded to the input/ folder if not already present.
forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044, :Tair => u"°C");
# get preset soil and initial temperature profile for Samoylov
soilprofile, tempprofile = CryoGrid.Presets.SamoylovDefault
initT = initializer(:T, tempprofile)
# choose grid with 5cm spacing
grid = CryoGrid.Presets.DefaultGrid_5cm
# build Tile from the given soil and temperature profiles
tile = CryoGrid.Presets.SoilHeatTile(TemperatureBC(forcings.Tair), GeothermalHeatFlux(0.053u"W/m^2"), soilprofile, initT, grid=grid)
# define time span (1 year)
tspan = (DateTime(2010,11,30),DateTime(2011,11,30))
u0, du0 = initialcondition!(tile, tspan)
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(tile, u0, tspan, savevars=(:T,))
# solve discretized system, saving every 3 hours;
out = @time solve(prob, saveat=3*3600.0, progress=true) |> CryoGridOutput;
zs = [2,7,12,22,32,42,50,100,500]u"cm"
cg = Plots.cgrad(:copper,rev=true)
plot(out.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false)
```
![Ts_output_freew](res/Ts_H_tair_freeW_2010-2011.png)
