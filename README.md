# CryoGrid.jl

Julia implementation of the CryoGrid land surface model using `DifferentialEquations.jl` and the [SciML](https://github.com/SciML)
package ecosystem.

Part of the broader research project: [Quantifying and explaining uncertainty in permafrost modeling under a warming climate](https://drive.google.com/file/d/1wB_EXtlO_PMXFSzZ-bRV8cg0a0DGDtAB/view?usp=sharing)

Author: Brian Groenke (brian.groenke@awi.de)

### Quick start

```julia
using CryoGrid
using CryoGrid.Models
using Dates
using Plots

forcings = loadforcings("input/FORCING_JSONfiles/FORCING_ULC_126_72.json", :Tair => u"°C");
# use air temperature as upper boundary forcing
tair = TimeSeriesForcing(ustrip.(u"K", forcings.data.Tair), forcings.timestamps, :Tair);
# basic 1-layer heat conduction model (defaults to free water freezing scheme)
model = Models.SoilHeat(TemperatureGradient(tair), SamoylovDefault)
# define time span
tspan = (DateTime(2010,1,1),DateTime(2010,12,31))
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(model,tspan)
# solve discretized system, saving every 6 hours;
# ROS3P is a third order Rosenbrock solver that should work well without a freeze curve.
out = @time solve(prob, ROS3P(), abstol=1e-2, saveat=6*3600.0) |> CryoGridOutput;
zs = [1:10...,20:10:100...]
cg = Plots.cgrad(:berlin,rev=true)
plot(out.soil.T[Z(zs)], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, dpi=150)

# Alternatively, we can use a van Genuchten freeze curve
model = Models.SoilHeat(TemperatureGradient(tair), SamoylovDefault, freezecurve=SFCC(VanGenuchten()))
# Set-up parameters
p = copy(model.pproto)
p.soil.α .= 4.0
p.soil.n .= 2.0
p.soil.Tₘ .= 273.15
tspan = (DateTime(2010,1,1),DateTime(2010,12,31))
prob = CryoGridProblem(model,tspan,p)
# stiff solvers don't work well with van Genuchten due to the ill-conditioned Jacobian;
# Thus, we use forward Euler instead
out = @time solve(prob, Euler(), dt=2*60.0, saveat=6*3600.0) |> CryoGridOutput;
zs = [1:10...,20:10:100...]
cg = Plots.cgrad(:berlin,rev=true)
plot(out.soil.T[Z(zs)], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, dpi=150)
```

Note that `SoilHeat` uses energy as the state variable by default. To use temperature as the state variable instead:

```julia
# Note that this will work with any freeze curve, here we use Westermann (2011).
# u"K" is the unit (Kelvin) for temperature, u"J" (Joules) represents energy.
# This is used in the specification of the Heat process type.
model = Models.SoilHeat(u"K", TemperatureGradient(tair), SamoylovDefault, freezecurve=SFCC(Westermann()))
```
