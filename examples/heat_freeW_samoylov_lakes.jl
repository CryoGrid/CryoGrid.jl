using CryoGrid
using Plots

forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044, :Tair => u"°C");
# define time span
tspan = (DateTime(2010,10,30),DateTime(2011,10,30))
# use air temperature as upper boundary forcing;
tair = TimeSeriesForcing(forcings.data.Tair, forcings.timestamps, :Tair);
T0 = values(tair[tspan[1]])[1]
tempprofile = TemperatureProfile(
    0.0u"m" => T0,
    1.0u"m" => -8.0u"°C",
    20.0u"m" => -10u"°C",
    1000.0u"m" => 1.0u"°C"
)
# soil profile: depth => (excess ice, natural porosity, saturation, organic fraction)
soilprofile = SoilProfile(
    0.0u"m" => HomogeneousMixture(por=0.80,sat=1.0,org=0.75), 
    0.1u"m" => HomogeneousMixture(por=0.80,sat=1.0,org=0.25),
    0.4u"m" => HomogeneousMixture(por=0.55,sat=1.0,org=0.25),
    3.0u"m" => HomogeneousMixture(por=0.50,sat=1.0,org=0.0),
    10.0u"m" => HomogeneousMixture(por=0.30,sat=1.0,org=0.0),
    100.0u"m" => HomogeneousMixture(por=0.10,sat=0.1,org=0.0),
);
initT = initializer(:T, tempprofile)
# Heat diffusion with temperature as prognostic state variable
# Construct soil layers
soil_layers = map(enumerate(soilprofile)) do (i, (depth, para))
    name = Symbol(:soil, i)
    depth => name => Soil(HeatBalance(), para=para)
end
# Build stratigraphy:
# @Stratigraphy macro lets us list multiple subsurface layers
strat = @Stratigraphy(
    -2.0*u"m" => Top(TemperatureGradient(tair)),
    -2.0u"m" => :lake => Lake(HeatBalance()),
    soil_layers...,
    998.0u"m" => Bottom(GeothermalHeatFlux(0.053u"J/s/m^2")),
);
# shift grid up by 2 m
modelgrid = Grid(CryoGrid.Presets.DefaultGrid_5cm .- 2.0u"m")
tile = Tile(strat, modelgrid, initT)
u0, du0 = initialcondition!(tile, tspan)
# construct CryoGridProblem with tile, initial condition, and timespan;
# we disable the default timestep limiter since we will use an adaptive solver.
prob = CryoGridProblem(tile, u0, tspan, savevars=(:T,))
# test step function (good to do debugging here)
tile(du0, u0, prob.p, prob.tspan[1])
@info "Running model"
out = @time solve(prob, Euler(), dt=300.0, saveat=24*3600.0, progress=true) |> CryoGridOutput;
# Plot it!
zs = [1,10,20,30,50,100,200,500,1000]u"cm"
cg = Plots.cgrad(:copper,rev=true);
plot(out.H[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Enthalpy", leg=false, dpi=150)
plot(out.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)
