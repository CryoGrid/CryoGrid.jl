# # [Soil heat with bucket water scheme](@id example5)
# In this example, we construct a `Tile` consisting of a soil column with (i) heat conduction
# and (ii) a bucket hydrology scheme. The rainfall data comes
# from the ERA5-Interim reanalysis product.

# Frist, load forcings and define boundary conditions.
using CryoGrid
forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044);
_, tempprofile = CryoGrid.Presets.SamoylovDefault;
initT = initializer(:T, tempprofile)
initsat = initializer(:sat, (l,state) -> state.sat .= l.para.sat);
upperbc = WaterHeatBC(SurfaceWaterBalance(forcings), TemperatureGradient(forcings.Tair, NFactor(0.5,0.9)));

# The `@Stratigraphy` macro is just a small convenience that automatically wraps the three subsurface layers in a tuple.
# It would be equivalent to use the `Stratigraphy` constructor directly and wrap these layers in a tuple or list.
strat = @Stratigraphy(
    0.0u"m" => Top(upperbc),
    0.0u"m" => Ground(MineralOrganic(por=0.80,sat=0.5,org=0.75), heat=HeatBalance(), water=WaterBalance(BucketScheme())),
    0.2u"m" => :middle => Ground(MineralOrganic(por=0.40,sat=0.75,org=0.10), heat=HeatBalance(), water=WaterBalance(BucketScheme())),
    2.0u"m" => Ground(MineralOrganic(por=0.10,sat=1.0,org=0.0), heat=HeatBalance(), water=WaterBalance(BucketScheme())),
    1000.0u"m" => Bottom(GeothermalHeatFlux(0.053u"W/m^2")),
);
modelgrid = CryoGrid.Presets.DefaultGrid_2cm
tile = Tile(strat, modelgrid, initT, initsat);

# Now we set up the problem and solve using the integrator interface.
tspan = (DateTime(2010,5,30),DateTime(2012,8,30))
u0, du0 = initialcondition!(tile, tspan)
prob = CryoGridProblem(tile, u0, tspan, savevars=(:T,:jw,:θw,:θwi), saveat=3*3600.0)
integrator = init(prob, CGEuler())
## Take a single step:
step!(integrator)
## ...then iterate over the remaining steps.
@time for i in integrator
    ## can add code here if necessary
end
out = CryoGridOutput(integrator.sol)

# Now let's check mass conservation for water.
water_added = values(sum(forcings.rainfall.(tspan[1]:Hour(3):tspan[2]).*u"m/s".*(3*3600.0u"s")))[1]
water_mass = Diagnostics.integrate(out.θwi, tile.grid)
Δwater = water_mass[end] - water_mass[1]

# Plot the results:
import Plots
zs = [1,5,9,15,21,25,33,55,65,75,100]u"cm"
cg = Plots.cgrad(:copper,rev=true);

# Temperature:
Plots.plot(out.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)

# Liquid water:
Plots.plot(out.θw[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Unfrozen water content", leg=false, size=(800,500), dpi=150)

# Saturation:
Plots.plot(out.sat[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Saturation", leg=false, size=(800,500), dpi=150)

# Runoff:
Plots.plot(out.top.runoff[1,:], ylabel="Runoff")