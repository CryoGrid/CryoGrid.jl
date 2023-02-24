using CryoGrid
using FreezeCurves
using FreezeCurves.Solvers
using Dates
using Plots

forcings = loadforcings(
    CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044,
    :Tair => u"°C",
    :rainfall => u"mm",
);
# use air temperature as upper boundary forcing;
tair = TimeSeriesForcing(forcings.data.Tair, forcings.timestamps, :Tair);
pr = TimeSeriesForcing(uconvert.(u"m/s", forcings.data.rainfall./3u"hr"), forcings.timestamps, :rainfall)
# define time span for simulation
tspan = (DateTime(2011,1,1),DateTime(2012,1,1))
T0 = values(tair[tspan[1]])[1]
tempprofile = TemperatureProfile(
    0.0u"m" => T0,
    1.0u"m" => -8.0u"°C",
    20.0u"m" => -10u"°C",
    1000.0u"m" => 1.0u"°C"
)
# soil profile: depth => (excess ice, natural porosity, saturation, organic fraction)
soilprofile = SoilProfile(
    0.0u"m" => HomogeneousMixture(por=0.80,sat=0.8,org=0.75), 
    0.1u"m" => HomogeneousMixture(por=0.80,sat=0.9,org=0.25),
    0.4u"m" => HomogeneousMixture(por=0.55,sat=0.95,org=0.25),
    3.0u"m" => HomogeneousMixture(por=0.50,sat=1.0,org=0.0),
    10.0u"m" => HomogeneousMixture(por=0.30,sat=1.0,org=0.0),
    100.0u"m" => HomogeneousMixture(por=0.10,sat=0.1,org=0.0),
);
initT = initializer(:T, tempprofile)
# initialize saturation to match soil profile
initsat = initializer(:sat, (l,p,state) -> state.sat .= l.para.sat)
# soil water retention curve and freeze curve
swrc = VanGenuchten(α=0.1, n=1.8)
sfcc = PainterKarra(ω=0.2, swrc=swrc)
# water flow: bucket scheme vs richard's eq
# waterflow = BucketScheme()
waterflow = RichardsEq(swrc=swrc)
# Heat diffusion with temperature as prognostic state variable
op = Heat.Diffusion(:H)
# @Stratigraphy macro lets us list multiple subsurface layers
strat = @Stratigraphy(
    -2.0*u"m" => Top(TemperatureGradient(tair), Rainfall(pr)),
    map(enumerate(soilprofile)) do (i, (depth, para))
        name = Symbol(:soil, i)
        depth => name => Soil(Coupled(WaterBalance(waterflow), HeatBalance(op, freezecurve=sfcc, dtlim=CryoGrid.MaxDelta(50u"kJ"))), para=para)
    end,
    1000.0u"m" => Bottom(GeothermalHeatFlux(0.053u"J/s/m^2")),
);
grid = CryoGrid.Presets.DefaultGrid_5cm
tile = Tile(strat, grid, initT, initsat);
u0, du0 = initialcondition!(tile, tspan)
prob = CryoGridProblem(tile, u0, tspan, saveat=3*3600, savevars=(:T,:θw,:θwi))
# note that this is currently quite slow since the Euler integrator takes very small time steps during the thawed season;
# expect it to take about 10-15 minutes per year on a typical workstation/laptop; minor speed-ups might be possible by tweaking the dt limiters
integrator = init(prob, Euler(), dt=60.0, saveat=3*3600.0)
@time while integrator.t < prob.tspan[end]
    # run the integrator forward in 10-day increments
    step!(integrator, 10*24*3600.0)
    t = convert_t(integrator.t)
    @info "t=$t, current dt=$(integrator.dt*u"s")"
end
out = CryoGridOutput(integrator.sol)
#out = @time solve(prob, Euler(), dt=60.0, saveat=3*3600.0, progress=true) |> CryoGridOutput;
# check mass conservation
water_added = values(sum(pr[tspan[1]:Hour(3):tspan[2]].tarray.*(3*3600.0u"s")))[1]
water_mass = Diagnostics.integrate(out.θwi, tile.grid)
Δwater = water_mass[end] - water_mass[1]
# Plot it!
zs = [1,5,10,15,20,30,40,50,100,150,200]u"cm"
cg = Plots.cgrad(:copper,rev=true);
plot(out.H[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Enthalpy", leg=false, dpi=150)
plot(out.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)
plot(out.sat[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Saturation", leg=false, size=(800,500), dpi=150)
