# CryoGrid.jl

Julia implementation of the CryoGrid land surface model using `DifferentialEquations.jl` and the [SciML](https://github.com/SciML)
package ecosystem.

Part of the broader research project: [Quantifying and explaining uncertainty in permafrost modeling under a warming climate](https://drive.google.com/file/d/1wB_EXtlO_PMXFSzZ-bRV8cg0a0DGDtAB/view?usp=sharing)

Author: Brian Groenke (brian.groenke@awi.de)

### Quick start

Single layer heat conduction model with free water freeze curve and air temperature upper boundary condition:

```julia
# load provided forcing data from Samoylov;
# The forcing file will be automatically downloaded to the input/ folder if not already present.
forcings = loadforcings(Models.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044, :Tair => u"°C");
# use air temperature (converted to Kelvin) as upper boundary forcing
tair = TimeSeriesForcing(ustrip.(u"K", forcings.data.Tair), forcings.timestamps, :Tair);
# basic 1-layer heat conduction model (defaults to free water freezing scheme)
model = Models.SoilHeat(TemperatureGradient(tair), SamoylovDefault)
# define time span (1 year)
tspan = (DateTime(2010,10,30),DateTime(2011,10,30))
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(model,tspan)
# solve discretized system, saving every 6 hours;
# Trapezoid on a discretized PDE is analogous to the well known Crank-Nicolson method.
out = @time solve(prob, Trapezoid(), abstol=1e-3, reltol=1e-4, saveat=6*3600.0, progress=true) |> CryoGridOutput;
zs = [1.0,5,10,20,30,50,100,500,1000]u"cm"
cg = Plots.cgrad(:copper,rev=true)
plot(out.soil.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, dpi=150)
```
![Ts_output_freew](res/Ts_H_tair_freeW_2010-2011.png)

Alternatively, we can use a van Genuchten freeze curve:

```julia
model = Models.SoilHeat(TemperatureGradient(tair), SamoylovDefault, freezecurve=SFCC(VanGenuchten()))
# Set-up parameters
p = copy(model.pproto)
p.soil.α .= 4.0
p.soil.n .= 2.0
p.soil.Tₘ .= 273.15
tspan = (DateTime(2010,10,30),DateTime(2011,10,30))
prob = CryoGridProblem(model,tspan,p)
# stiff solvers don't work well with van Genuchten due to the ill-conditioned Jacobian;
# We can just forward Euler instead.
out = @time solve(prob, Euler(), dt=120.0, saveat=6*3600.0, progress=true) |> CryoGridOutput;
plot(out.soil.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, dpi=150)
```
![Ts_output_vgfc](res/Ts_H_tair_vg_2010-2011.png)

Note that `SoilHeat` uses energy as the state variable by default. To use temperature as the state variable instead:

```julia
# Note that this will work with any freeze curve, here we use Westermann (2011).
# :T is the variable name for temperature, :H represents enthalpy/energy.
# This is used in the specification of the Heat process type.
model = Models.SoilHeat(:T, TemperatureGradient(tair), SamoylovDefault, freezecurve=SFCC(Westermann()))
```
