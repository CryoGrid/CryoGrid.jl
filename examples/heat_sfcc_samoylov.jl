using CryoGrid

forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044);
# "simple" heat conduction model w/ 5 cm grid spacing
grid = CryoGrid.Presets.DefaultGrid_5cm
soilprofile, tempprofile = CryoGrid.Presets.SamoylovDefault
initT = initializer(:T, tempprofile)
sfcc = PainterKarra(swrc=VanGenuchten(α=0.1, n=2.0))
upperbc = TemperatureGradient(forcings.Tair)
lowerbc = GeothermalHeatFlux(0.053u"W/m^2")
tile = CryoGrid.Presets.SoilHeatTile(:H, upperbc, lowerbc, soilprofile, initT; grid=grid, freezecurve=sfcc)
# define time span
tspan = (DateTime(2010,11,30),DateTime(2011,11,30))
u0, du0 = initialcondition!(tile, tspan)
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(tile, u0, tspan, saveat=3*3600.0, savevars=(:T,))
@info "Running model"
out = @time solve(prob, SSPRK22(), dt=300.0, saveat=3*3600.0, progress=true) |> CryoGridOutput;

# Plot it!
import CairoMakie

zs = [5,10,15,20,25,30,40,50,100,500]u"cm"
Diagnostics.plot_at_depths(:T, out, zs, ylabel="Temperature (°C)")
