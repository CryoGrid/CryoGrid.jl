using CryoGrid
using Plots

# Custom grid (though actually pretty much the same as CryoGridModels.DefaultGrid_2cm)
const gridvals = vcat([0:0.02:2...,2.05:0.05:4.0...,
	4.1:0.1:10...,10.2:0.2:20...,21:1:30...,
	35:5:50...,60:10:100...,200:100:1000...]...)
# soil profile: depth => (excess ice, natural porosity, saturation, organic fraction)
soilprofile = SoilProfile(
	0.0u"m" => SoilProperties(χ=0.0,ϕ=0.80,θ=1.0,ω=0.75), #(θw=0.80,θm=0.05,θo=0.15,ϕ=0.80),
	0.1u"m" => SoilProperties(χ=0.0,ϕ=0.80,θ=1.0,ω=0.25), #(θw=0.80,θm=0.15,θo=0.05,ϕ=0.80),
	0.4u"m" => SoilProperties(χ=0.30,ϕ=0.55,θ=1.0,ω=0.25), #(θw=0.80,θm=0.15,θo=0.05,ϕ=0.55),
	3.0u"m" => SoilProperties(χ=0.0,ϕ=0.50,θ=1.0,ω=0.0), #(θw=0.50,θm=0.50,θo=0.0,ϕ=0.50),
	10.0u"m" => SoilProperties(χ=0.0,ϕ=0.30,θ=1.0,ω=0.0), #(θw=0.30,θm=0.70,θo=0.0,ϕ=0.30),
)
# mid-winter temperature profile
tempprofile = TempProfile(
    0.01u"m" => -15.9u"°C",
    0.05u"m" => -15.49u"°C",
    0.10u"m" => -15.32u"°C",
    0.20u"m" => -14.44u"°C",
    0.30u"m" => -14.18u"°C",
    0.40u"m" => -13.50u"°C",
    1000.0u"m" => 10.2u"°C",
)
forcings = loadforcings(CryoGridModels.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044, :Tair => u"°C");
# convert Tair to Kelvin
Tair = TimeSeriesForcing(ustrip.(u"K", forcings.data.Tair), forcings.timestamps, :Tair);
# assume other forcings don't (yet) have units
pr   = TimeSeriesForcing(forcings.data.pressure, forcings.timestamps, :p);
q    = TimeSeriesForcing(forcings.data.q, forcings.timestamps, :q);
wind = TimeSeriesForcing(forcings.data.wind, forcings.timestamps, :wind);
Lin  = TimeSeriesForcing(forcings.data.Lin, forcings.timestamps, :Lin);
Sin  = TimeSeriesForcing(forcings.data.Sin, forcings.timestamps, :Sin);
z = 2.;    # height [m] for which the forcing variables (Temp, humidity, wind, pressure) are provided
tspan = (DateTime(2010,1,1), DateTime(2011,1,1))
strat = Stratigraphy(
    -2.0u"m" => Top(SurfaceEnergyBalance(Tair,pr,q,wind,Lin,Sin,z)),
    # Note: You can change Heat{:H} to Heat{:T} to use temperature as the prognostic state variable.
    0.0u"m" => Ground(:soil, Soil(soilprofile), Heat{:H}(tempprofile, freezecurve=SFCC(VanGenuchten()))),
    1000.0u"m" => Bottom(GeothermalHeatFlux(0.053u"J/s/m^2"))
);
grid = Grid(gridvals);
model = CryoGridSetup(strat,grid);
# define time span
tspan = (DateTime(2010,10,30),DateTime(2011,10,30))
p = copy(model.pproto)
p.soil.α .= 4.0
p.soil.n .= 2.0
p.soil.Tₘ .= 273.15 # K
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(model,tspan,p)
# OPTIONAL: solve for short time period to 
# solve with forward Euler (w/ CFL) and construct CryoGridOutput from solution
out = @time solve(prob, Euler(), dt=2*60.0, callback=CFLStepLimiter(model), saveat=24*3600.0, progress=true) |> CryoGridOutput;
# Plot it!
zs = [1:10...,20:10:100...]
cg = Plots.cgrad(:copper,rev=true);
plot(out.soil.H[Z(zs)], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Enthalpy", leg=false, dpi=150)
plot(out.soil.T[Z(zs)] .- 273.15, color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)
