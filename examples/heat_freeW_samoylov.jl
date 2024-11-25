# # [Soil heat with free water freeze curve](@id example1)
# In this example, we use the preset `SoilHeatTile` to construct
# a `Tile` consisting of a soil column with heat conduction
# forced using air temperatures from Samoylov Island.
# The enthalpy-based `HeatBalance` process defaults to the so-called
# "free water" freezing characteristic which assumes that water only
# freezes and thaws at a melting temperature of 0°C.

# First we load the built-in forcing file from Nitzbon et al. 2020 (CryoGrid 3). Note that this will
# download the forcing files from the AWI NextCloud if they are not already present in the `input/` folder.
using CryoGrid
forcings = loadforcings(CryoGrid.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044);

# We use a simple 5-layer stratigraphy suitable for Samoylov. This is based on the
# profile provided in `Presets` but uses the default "free water" freezing characteristic
# defined on `SimpleSoil`.
soilprofile = SoilProfile(
    0.0u"m" => SimpleSoil(; por=0.80, org=0.75),
    0.1u"m" => SimpleSoil(; por=0.80, org=0.25),
    0.4u"m" => SimpleSoil(; por=0.55, org=0.25),
    3.0u"m" => SimpleSoil(; por=0.50, org=0.0),
    10.0u"m" => SimpleSoil(; por=0.30, org=0.0),
);

# We construct a state variable initializer for temperature `T` from the temperature profile preset for Samoylov.
initT = initializer(:T, CryoGrid.Presets.SamoylovDefault.tempprofile);

# We choose the default grid discretization with 5 cm spacing at the surface.
grid = CryoGrid.DefaultGrid_5cm;

# Now we construct the Tile using the built-in model configuration `SoilHeatTile` which defines a
# standalone soil straigraphy with only heat conduction and no water flow.
tile = CryoGrid.Presets.SoilHeatTile(
    :H,
    TemperatureBC(Input(:Tair), NFactor(nf=Param(0.6), nt=Param(0.9))),
    GeothermalHeatFlux(0.053u"W/m^2"),
    soilprofile,
    forcings,
    initT;
    grid=grid
);

# Here we define the time span:
tspan = (DateTime(2010,12,31),DateTime(2011,12,31));

# Evaluate the initial condition
u0, du0 = initialcondition!(tile, tspan);

# Here we construct a CryoGridProblem with tile, initial condition, and timespan.
prob = CryoGridProblem(tile, u0, tspan, saveat=24*3600.0, savevars=(:T,));

# Solve the configured problem with the built-in forward Euler method.
# note that, due to compile time, this may take 1-2 minutes when executed in a fresh Julia
# session. Subsequent solves will be much faster.
sol = @time solve(prob);
out = CryoGridOutput(sol)

# Now we plot the reuslts!
import Plots
zs = [1,10,20,30,50,100,200,500,1000]u"cm"
cg = Plots.cgrad(:copper,rev=true);
Plots.plot(out.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)

using BenchmarkTools
@profview @btime $tile($du0, $u0, $prob.p, $prob.tspan[1])
