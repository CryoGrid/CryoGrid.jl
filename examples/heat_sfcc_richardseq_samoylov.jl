using CryoGrid

forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044);
# define time span for simulation
tspan = (DateTime(2011,1,1),DateTime(2012,1,1))
T0 = values(forcings.Tair(tspan[1]))[1]
tempprofile = TemperatureProfile(
    0.0u"m" => T0*u"°C",
    1.0u"m" => -8.0u"°C",
    20.0u"m" => -10u"°C",
    1000.0u"m" => 1.0u"°C"
)
initT = initializer(:T, tempprofile)
# initialize saturation to match soil profile
initsat = initializer(:sat, (l,state) -> state.sat .= l.para.sat)
# soil water retention curve and freeze curve
swrc = VanGenuchten(α=0.1, n=1.8)
sfcc = PainterKarra(ω=0.0, swrc=swrc)
# water flow: bucket scheme vs richard's eq
# waterflow = BucketScheme()
waterflow = RichardsEq(swrc=swrc)
# Enthalpy-based heat diffusion with high accuracy Newton-based solver for inverse enthalpy mapping
heatop = Heat.EnthalpyForm(SFCCNewtonSolver())
upperbc = WaterHeatBC(SurfaceWaterBalance(rainfall=forcings.rainfall), TemperatureGradient(forcings.Tair, NFactor()))
# We will use a simple stratigraphy with 3 subsurface soil layers
# Note that the @Stratigraphy macro lets us list multiple subsurface layers
strat = @Stratigraphy(
    -2.0u"m" => Top(upperbc),
    0.0u"m" => :topsoil => HomogeneousSoil(MineralOrganic(por=0.80,sat=0.7,org=0.75), heat=HeatBalance(heatop, freezecurve=sfcc), water=WaterBalance(RichardsEq(;swrc))),
    0.2u"m" => :subsoil => HomogeneousSoil(MineralOrganic(por=0.40,sat=0.8,org=0.10), heat=HeatBalance(heatop, freezecurve=sfcc), water=WaterBalance(RichardsEq(;swrc))),
    2.0u"m" => :substrat => HomogeneousSoil(MineralOrganic(por=0.10,sat=1.0,org=0.0), heat=HeatBalance(heatop, freezecurve=sfcc), water=WaterBalance(RichardsEq(;swrc))),
    1000.0u"m" => Bottom(GeothermalHeatFlux(0.053u"W/m^2"))
);
grid = CryoGrid.Presets.DefaultGrid_2cm
tile = Tile(strat, grid, initT, initsat);
u0, du0 = initialcondition!(tile, tspan)
prob = CryoGridProblem(tile, u0, tspan, saveat=3*3600, savevars=(:T,:θw,:θwi,:kw))
# note that this is currently somewhat slow since the integrator must take very small time steps during the thawed season;
# expect it to take about 3-5 minutes per year on a typical workstation/laptop;
# minor speed-ups might be possible by tweaking the dt limiters or by using the SFCCPreSolver
integrator = init(prob, Euler(), dt=1.0, saveat=3*3600.0)
# take one step to check if it's working
step!(integrator)
# we can use the `getstate` function to construct the current Tile state from the integrator
state = getstate(integrator)
# check that all water fluxes are near zero since we're starting in frozen conditions
@assert all(isapprox.(0.0, state.topsoil.jw, atol=1e-14))
# run the integrator forward in time until the end of the tspan
@time while integrator.t < prob.tspan[end]
    @assert all(isfinite.(integrator.u))
    @assert all(0 .<= integrator.u.sat .<= 1)
    # run the integrator forward in daily increments
    step!(integrator, 24*3600.0)
    t = convert_t(integrator.t)
    @info "t=$t, current dt=$(integrator.dt*u"s")"
end

out = CryoGridOutput(integrator.sol)

# check mass conservation; TODO: need to track surface water runoff to close the water balance
water_added = values(sum(upreferred.(forcings.rainfall.(tspan[1]:Hour(3):tspan[2]).*u"m/s".*3u"hr")))[1]
water_mass = Diagnostics.integrate(out.θwi, tile.grid)
water_resid = water_mass[end] - water_mass[1] - water_added*1u"m^2"

# Plot it!
import Plots

zs = [1,5,10,15,20,30,40,50,100,150,200]u"cm"
cg = Plots.cgrad(:copper,rev=true);
Plots.plot(out.H[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Enthalpy", leg=false, dpi=150)
Plots.plot(out.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)
Plots.plot(out.sat[Z(1:10)], color=cg[LinRange(0.0,1.0,10)]', ylabel="Saturation", leg=false, size=(800,500), dpi=150)
