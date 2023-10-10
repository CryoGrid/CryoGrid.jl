# # [Soil heat with SFCC](@id example3)
# In this example, we use the preset `SoilHeatTile` to construct
# a `Tile` consisting of a soil column with heat conduction
# forced using air temperatures from Samoylov Island. We use the
# SFCC formulation of Painter and Karra (2014). For the purpose
# of demonstration, we use the apparent heat capacity form of the
# heat equation in this example (i.e. [`Heat.MOLTemperature`](@ref)).

using CryoGrid

# First we set up the model:
forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA5_fitted_daily_1979_2020);
grid = CryoGrid.Presets.DefaultGrid_5cm
soilprofile, tempprofile = CryoGrid.Presets.SamoylovDefault
initT = initializer(:T, tempprofile)
sfcc = PainterKarra(swrc=VanGenuchten(α=0.1, n=2.0))
upperbc = TemperatureGradient(forcings.Tair, NFactor(nf=0.6))
lowerbc = GeothermalHeatFlux(0.053u"W/m^2")
tile = CryoGrid.Presets.SoilHeatTile(upperbc, lowerbc, soilprofile, initT; grid=grid, freezecurve=sfcc)
tspan = (DateTime(2010,10,30),DateTime(2011,10,30))
u0, du0 = initialcondition!(tile, tspan)
prob = CryoGridProblem(tile, u0, tspan, saveat=3*3600.0, savevars=(:T,));

# ... then solve it with the built-in forward Euler integrator.
sol = @time solve(prob);
out = CryoGridOutput(sol)

# Finally, plot the resulting temperatures.
import Plots
zs = [5,10,15,20,25,30,40,50,100,500]u"cm"
Diagnostics.plot_at_depths(:T, out, zs, ylabel="Temperature (°C)")
