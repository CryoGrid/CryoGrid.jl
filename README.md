# CryoGrid.jl

[![][docs-dev-img]][docs-dev-url]

[docs-dev-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-dev-url]: https://cryogrid.github.io/CryoGrid.jl/dev/

Julia implementation of the CryoGrid permafrost model using [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) and the [SciML](https://github.com/SciML)
package ecosystem.

Part of the broader research project: [Quantifying and explaining uncertainty in permafrost modeling under a warming climate](https://drive.google.com/file/d/1wB_EXtlO_PMXFSzZ-bRV8cg0a0DGDtAB/view?usp=sharing)

Author: Brian Groenke (brian.groenke@awi.de)

### Installation

`CryoGrid.jl` can be installed via the Julia package manager:

```
add CryoGrid
```

or equivalently in code/REPL:

```julia
import Pkg
Pkg.add("CryoGrid")
```

### Quick start

Single layer heat conduction model with free water freeze curve and simple air temperature upper boundary condition:

```julia
using CryoGrid
using Plots: plot

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
