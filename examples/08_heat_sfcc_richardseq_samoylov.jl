# # [Coupled soil heat and water transport](@id example8)
# In this example, we construct a `Tile` from a Stratigraphy of three soil layers
# with coupled water and heat transport. We use the Richards-Richardson equation for
# unsaturated flow in porous media to account for flow due to capillary suction and
# pressure gradients.

using CryoGrid
forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044);
tspan = (DateTime(2011,1,1),DateTime(2011,12,31))
T0 = values(forcings.Tair(tspan[1]))[1]
tempprofile = TemperatureProfile(
    0.0u"m" => T0*u"°C",
    1.0u"m" => -8.0u"°C",
    20.0u"m" => -10u"°C",
    1000.0u"m" => 1.0u"°C"
)
initT = initializer(:T, tempprofile);

# Here we conigure the water retention curve and freeze curve. The van Genuchten parameters coorespond to that
# which would be reasonable for a silty soil.
swrc = VanGenuchten(α=0.1, n=1.8)
sfcc = PainterKarra(ω=0.0, swrc=swrc)
waterflow = RichardsEq(swrc=swrc);

# We use the enthalpy-based heat diffusion with high accuracy Newton-based solver for inverse enthalpy mapping
heatop = Heat.EnthalpyForm(SFCCNewtonSolver())
upperbc = WaterHeatBC(SurfaceWaterBalance(rainfall=forcings.rainfall), TemperatureGradient(forcings.Tair, NFactor(nf=0.6, nt=0.9)));

# We will use a simple stratigraphy with three subsurface soil layers.
# Note that the @Stratigraphy macro lets us list multiple subsurface layers without wrapping them in a tuple.
strat = @Stratigraphy(
    -2.0u"m" => Top(upperbc),
    0.0u"m" => :topsoil => Ground(MineralOrganic(por=0.80,sat=0.7,org=0.75), heat=HeatBalance(heatop, freezecurve=sfcc), water=WaterBalance(RichardsEq(;swrc))),
    0.2u"m" => :subsoil => Ground(MineralOrganic(por=0.40,sat=0.8,org=0.10), heat=HeatBalance(heatop, freezecurve=sfcc), water=WaterBalance(RichardsEq(;swrc))),
    2.0u"m" => :substrat => Ground(MineralOrganic(por=0.10,sat=1.0,org=0.0), heat=HeatBalance(heatop, freezecurve=sfcc), water=WaterBalance(RichardsEq(;swrc))),
    1000.0u"m" => Bottom(GeothermalHeatFlux(0.053u"W/m^2"))
);
grid = CryoGrid.Presets.DefaultGrid_2cm
tile = Tile(strat, grid, initT);
u0, du0 = initialcondition!(tile, tspan)
prob = CryoGridProblem(tile, u0, tspan, saveat=3*3600, savevars=(:T,:θw,:θwi,:kw));

# This is currently somewhat slow since the integrator must take very small time steps during the thawed season;
# expect it to take about 3-5 minutes per year on a typical workstation/laptop.
# Minor speed-ups might be possible by tweaking the dt limiters or by using the `SFCCPreSolver` for the freeze curve.
integrator = init(prob, CGEuler())

# Here we take just one step to check if it's working.
step!(integrator)
@assert all(isfinite.(integrator.u))

# We can use the `getstate` function to construct the current Tile state from the integrator.
# We then check that all water fluxes are near zero since we're starting in frozen conditions.
state = getstate(integrator);
@assert all(isapprox.(0.0, state.topsoil.jw, atol=1e-14))

# Run the integrator forward in time until the end of the tspan:
@time while integrator.t < prob.tspan[end]
    @assert all(isfinite.(integrator.u))
    @assert all(0 .<= integrator.u.sat .<= 1)
    ## run the integrator forward in daily increments
    step!(integrator, 24*3600.0)
    t = convert_t(integrator.t)
    @info "t=$t, current dt=$(integrator.dt*u"s")"
end;
out = CryoGridOutput(integrator.sol)

# Check mass conservation...
water_added = values(sum(upreferred.(forcings.rainfall.(tspan[1]:Hour(3):tspan[2]).*u"m/s".*3u"hr")))[1]
water_mass = Diagnostics.integrate(out.θwi, tile.grid)
water_resid = water_mass[end] - water_mass[1] - water_added*1u"m^2"

# Plot it!
import Plots
zs = [1,5,10,15,20,30,40,50,100,150,200]u"cm"
cg = Plots.cgrad(:copper,rev=true);
p1 = Plots.plot(out.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)
p2 = Plots.plot(out.sat[Z(1:10)], color=cg[LinRange(0.0,1.0,10)]', ylabel="Saturation", leg=false, size=(800,500), dpi=150)
Plots.plot(p1, p2, size=(1200,400))
