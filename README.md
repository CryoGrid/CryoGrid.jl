# CryoGrid.jl

[![][docs-dev-img]][docs-dev-url]

[docs-dev-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-dev-url]: https://cryogrid.github.io/CryoGrid.jl/dev/

Julia implementation of the CryoGrid permafrost model using `DifferentialEquations.jl` and the [SciML](https://github.com/SciML)
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
Pkg.add(url="https://gitlab.awi.de/sparcs/cryogrid/cryogridjulia")
```

### Quick start

Single layer heat conduction model with free water freeze curve and air temperature upper boundary condition:

```julia
using CryoGrid
using Plots

# load provided forcing data from Samoylov;
# The forcing file will be automatically downloaded to the input/ folder if not already present.
forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044);
# get preset soil and initial temperature profile for Samoylov
soilprofile, tempprofile = CryoGrid.Presets.SamoylovDefault
initT = initializer(:T, tempprofile)
# choose grid with 5cm spacing
grid = CryoGrid.Presets.DefaultGrid_5cm
# basic 1-layer heat conduction model (defaults to free water freezing scheme)
tile = CryoGrid.Presets.SoilHeatTile(
    TemperatureGradient(forcings.Tair),
    GeothermalHeatFlux(0.053u"W/m^2"),
    soilprofile,
    initT,
    grid=grid
)
# define time span (1 year)
tspan = (DateTime(2010,11,30),DateTime(2011,11,30))
u0, du0 = initialcondition!(tile, tspan)
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(tile, u0, tspan, savevars=(:T,))
# solve discretized system, saving every 3 hours;
# Trapezoid on a discretized PDE is analogous to the well known Crank-Nicolson method.
out = @time solve(prob, Trapezoid(), saveat=3*3600.0, progress=true) |> CryoGridOutput;
zs = [2,7,12,22,32,42,50,100,500]u"cm"
cg = Plots.cgrad(:copper,rev=true)
plot(out.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false)
```
![Ts_output_freew](res/Ts_H_tair_freeW_2010-2011.png)

Alternatively, we can use a Dall'Amico freeze curve:

```julia
sfcc = DallAmico(swrc=VanGenuchten(α=0.02, n=1.8)) # silt/clay-like freeze curve
tile2 = CryoGrid.Presets.SoilHeatTile(
    TemperatureGradient(forcings.Tair),
    GeothermalHeatFlux(0.053u"W/m^2"),
    soilprofile,
    initT,
    grid=grid,
    freezecurve=sfcc
)
u0, du0 = initialcondition!(tile2, tspan)
# CryoGrid front-end for ODEProblem
prob2 = CryoGridProblem(tile2, u0, tspan, savevars=(:T,))
# stiff solvers don't work well with Dall'Amico due to the ill-conditioned Jacobian;
# We can just forward Euler instead.
out2 = @time solve(prob2, Euler(), dt=300.0, saveat=3*3600.0, progress=true) |> CryoGridOutput;
plot(out2.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false)
```

![Ts_output_freew](res/Ts_H_tair_vg_2010-2011.png)

Note that `SoilHeatTile` uses energy as the state variable by default. To use temperature as the state variable instead:

```julia
# :T is the variable name for temperature, :H represents enthalpy/energy.
# This is used in the specification of the HeatBalance process type.
# While this will work with any freeze curve, here we use Westermann (2011) as an example.
model = CryoGrid.Presets.SoilHeatTile(:T, TemperatureGradient(forcings.Tair), soilprofile, freezecurve=SFCC(Westermann()))
```
