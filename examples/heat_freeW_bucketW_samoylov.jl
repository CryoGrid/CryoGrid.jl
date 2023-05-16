using CryoGrid
using Dates
using Plots

forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044);
grid = CryoGrid.Presets.DefaultGrid_5cm
_, tempprofile = CryoGrid.Presets.SamoylovDefault
initT = initializer(:T, tempprofile)
initsat = initializer(:sat, (l,p,state) -> state.sat .= l.para.sat)
# We will use a simple stratigraphy with 3 subsurface soil layers
# Note that the @Stratigraphy macro lets us list multiple subsurface layers
strat = @Stratigraphy(
    -2.0u"m" => Top(TemperatureGradient(forcings.Tair), Rainfall(forcings.rainfall)),
    0.0u"m" => :topsoil => HomogeneousSoil(MineralOrganic(por=0.80,sat=0.7,org=0.75), heat=HeatBalance(), water=WaterBalance(BucketScheme())),
    0.2u"m" => :subsoil => HomogeneousSoil(MineralOrganic(por=0.40,sat=0.8,org=0.10), heat=HeatBalance(), water=WaterBalance(BucketScheme())),
    2.0u"m" => :substrat => HomogeneousSoil(MineralOrganic(por=0.10,sat=1.0,org=0.0), heat=HeatBalance(), water=WaterBalance(BucketScheme())),
    1000.0u"m" => Bottom(GeothermalHeatFlux(0.053u"W/m^2"))
);
tile = Tile(strat, grid, initT, initsat);
# define time span
tspan = (DateTime(2011,10,30),DateTime(2012,10,30))
u0, du0 = initialcondition!(tile, tspan)
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(tile, u0, tspan, savevars=(:T,:jw,:θw,:θwi), saveat=3*3600.0)
# solve with forward Euler and initial timestep of 5 minutes
integrator = init(prob, Euler(), dt=300.0, saveat=3*3600.0)
@time for i in integrator
    @assert all(isfinite.(integrator.u))
end
out = CryoGridOutput(integrator.sol)
# check mass conservation
water_added = values(sum(forcings.rainfall.(tspan[1]:Hour(3):tspan[2]).*u"m/s".*(3*3600.0u"s")))[1]
water_mass = Diagnostics.integrate(out.θwi, tile.grid)
Δwater = water_mass[end] - water_mass[1]
# Plot it!
zs = [1,5,10,15,20,30,40,50,100,150,200]u"cm"
cg = Plots.cgrad(:copper,rev=true);
plot(out.H[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Enthalpy", leg=false, dpi=150)
plot(out.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)
plot(out.θw[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Unfrozen water content", leg=false, size=(800,500), dpi=150)
plot(out.sat[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Saturation", leg=false, size=(800,500), dpi=150)
