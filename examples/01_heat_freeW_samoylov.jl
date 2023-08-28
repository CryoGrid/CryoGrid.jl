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
forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044);

# We use the provided default soil and temperature profiles for Samoylov.
soilprofile, tempprofile = CryoGrid.Presets.SamoylovDefault;

# We choose the default grid discretization with 5 cm spacing at the surface.
grid = CryoGrid.Presets.DefaultGrid_5cm;

# We construct a state variable initializer for temperature `T` from the temperature profile.
initT = initializer(:T, tempprofile)
tile = CryoGrid.Presets.SoilHeatTile(
    :H,
    TemperatureGradient(forcings.Tair),
    GeothermalHeatFlux(0.053u"W/m^2"),
    soilprofile,
    initT;
    grid=grid
);

# Here we define the time span:
tspan = (DateTime(2010,10,30),DateTime(2011,10,30));

# Evaluate the initial condition
u0, du0 = initialcondition!(tile, tspan);

# Here we construct a CryoGridProblem with tile, initial condition, and timespan;
# we disable the default timestep limiter since we will use an adaptive solver.
prob = CryoGridProblem(tile, u0, tspan, saveat=24*3600.0, savevars=(:T,:jH), step_limiter=nothing);

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
