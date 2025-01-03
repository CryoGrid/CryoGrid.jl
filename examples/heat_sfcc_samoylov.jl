# # [Soil heat with SFCC](@id example3)
# In this example, we use the preset `SoilHeatTile` to construct
# a `Tile` consisting of a soil column with heat conduction
# forced using n-factor scaled air temperatures from Samoylov Island.

using CryoGrid

# First we set up the soil heat model. Note that the default soil profile
# for Samoylov already includes the appropriate freeze curves.
forcings = loadforcings(CryoGrid.Forcings.Samoylov_ERA5_fitted_daily_1979_2020);
grid = CryoGrid.DefaultGrid_5cm
soilprofile, tempprofile = CryoGrid.SamoylovDefault
initT = initializer(:T, tempprofile)
upperbc = TemperatureBC(Input(:Tair), NFactor(nf=0.6, nt=0.9))
lowerbc = GeothermalHeatFlux(0.053u"W/m^2")
tile = CryoGrid.SoilHeatTile(upperbc, lowerbc, soilprofile, forcings, initT; grid=grid)
tspan = (DateTime(2010,10,30),DateTime(2011,10,30))
u0, du0 = @time initialcondition!(tile, tspan)
prob = CryoGridProblem(tile, u0, tspan, saveat=3*3600.0, savevars=(:T,));

# ... then solve it with the built-in forward Euler integrator.
sol = @time solve(prob);
out = CryoGridOutput(sol)

# Finally, plot the resulting temperatures.
import Plots
zs = [5,10,15,20,25,30,40,50,100,500]u"cm"
Diagnostics.plot_at_depths(:T, out, zs, ylabel="Temperature (°C)")
