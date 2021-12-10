# CryoGrid.jl

[![][docs-dev-img]][docs-dev-url]

[docs-dev-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-dev-url]: https://cryogrid.github.io/CryoGrid.jl/dev/

Julia implementation of the CryoGrid land surface model using `DifferentialEquations.jl` and the [SciML](https://github.com/SciML)
package ecosystem.

Part of the broader research project: [Quantifying and explaining uncertainty in permafrost modeling under a warming climate](https://drive.google.com/file/d/1wB_EXtlO_PMXFSzZ-bRV8cg0a0DGDtAB/view?usp=sharing)

Author: Brian Groenke (brian.groenke@awi.de)

### Installation

`CryoGrid.jl` can be installed via the Julia package manager:

```
add https://gitlab.awi.de/sparcs/cryogrid/cryogridjulia
```

or equivalently in code/REPL:

```julia
import Pkg
Pkg.add(["https://gitlab.awi.de/sparcs/cryogrid/cryogridjulia"])
```

### Quick start

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
model = CryoGrid.Presets.SoilHeatColumn(TemperatureGradient(tair), soilprofile)
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
model = CryoGrid.Presets.SoilHeatColumn(TemperatureGradient(tair), soilprofile, freezecurve=SFCC(DallAmico()))
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

Note that `SoilHeatColumn` uses energy as the state variable by default. To use temperature as the state variable instead:

```julia
# :T is the variable name for temperature, :H represents enthalpy/energy.
# This is used in the specification of the Heat process type.
# While this will work with any freeze curve, here we use Westermann (2011) as an example.
model = CryoGrid.Presets.SoilHeatColumn(:T, TemperatureGradient(tair), soilprofile, freezecurve=SFCC(Westermann()))
```
