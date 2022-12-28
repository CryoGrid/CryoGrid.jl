using CryoGrid
using Plots

# Custom grid (though actually pretty much the same as CryoGrid.Presets.DefaultGrid_2cm)
const gridvals = vcat([0:0.02:2...,2.05:0.05:4.0...,
	4.1:0.1:10...,10.2:0.2:20...,21:1:30...,
	35:5:50...,60:10:100...,200:100:1000...]...)u"m"
# soil profile: depth => (excess ice, natural porosity, saturation, organic fraction)
soilprofile = SoilProfile(
    0.0u"m" => HomogeneousMixture(por=0.80,sat=1.0,org=0.75), #(θwi=0.80,θm=0.05,θo=0.15,ϕ=0.80),
    0.1u"m" => HomogeneousMixture(por=0.80,sat=1.0,org=0.25), #(θwi=0.80,θm=0.15,θo=0.05,ϕ=0.80),
    0.4u"m" => HomogeneousMixture(por=0.55,sat=1.0,org=0.25), #(θwi=0.80,θm=0.15,θo=0.05,ϕ=0.55),
    3.0u"m" => HomogeneousMixture(por=0.50,sat=1.0,org=0.0), #(θwi=0.50,θm=0.50,θo=0.0,ϕ=0.50),
    10.0u"m" => HomogeneousMixture(por=0.30,sat=1.0,org=0.0), #(θwi=0.30,θm=0.70,θo=0.0,ϕ=0.30),
);
# mid-winter temperature profile
tempprofile = CryoGrid.Presets.SamoylovDefault.tempprofile
forcings = loadforcings(
    CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044,
    :Tair => u"°C",
    :pressure => u"kPa",
    :wind => u"m/s",
    :Lin => u"W/m^2",
    :Sin => u"W/m^2",
);
Tair = TimeSeriesForcing(forcings.data.Tair, forcings.timestamps, :Tair);
pr   = TimeSeriesForcing(forcings.data.pressure, forcings.timestamps, :pr);
q    = TimeSeriesForcing(forcings.data.q, forcings.timestamps, :q);
wind = TimeSeriesForcing(forcings.data.wind, forcings.timestamps, :wind);
Lin  = TimeSeriesForcing(forcings.data.Lin, forcings.timestamps, :Lin);
Sin  = TimeSeriesForcing(forcings.data.Sin, forcings.timestamps, :Sin);
z = 2.;    # height [m] for which the forcing variables (Temp, humidity, wind, pressure) are provided
tspan = (DateTime(2010,1,1), DateTime(2011,1,1))
soilprofile, tempprofile = CryoGrid.Presets.SamoylovDefault
initT = initializer(:T, tempprofile)
# @Stratigraphy macro lets us list multiple subsurface layers
strat = @Stratigraphy(
    -z*u"m" => Top(SurfaceEnergyBalance(Tair, pr, q,wind, Lin, Sin, z, solscheme=SEB.Iterative(), stabfun=SEB.HøgstrømSHEBA())),
    soilprofile[1].depth => :soil1 => Soil(HeatBalance(:H, freezecurve=SFCC(DallAmico())), para=soilprofile[1].value),
    soilprofile[2].depth => :soil2 => Soil(HeatBalance(:H, freezecurve=SFCC(DallAmico())), para=soilprofile[2].value),
    soilprofile[3].depth => :soil3 => Soil(HeatBalance(:H, freezecurve=SFCC(DallAmico())), para=soilprofile[3].value),
    soilprofile[4].depth => :soil4 => Soil(HeatBalance(:H, freezecurve=SFCC(DallAmico())), para=soilprofile[4].value),
    soilprofile[5].depth => :soil5 => Soil(HeatBalance(:H, freezecurve=SFCC(DallAmico())), para=soilprofile[5].value),
    1000.0u"m" => Bottom(GeothermalHeatFlux(0.053u"J/s/m^2")),
);
grid = Grid(gridvals);
tile = Tile(strat, grid, initT);
# define time span
tspan = (DateTime(2010,10,30),DateTime(2011,10,30))
u0, du0 = initialcondition!(tile, tspan)
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(tile, u0, tspan, savevars=(:T,))
# solve with forward Euler and construct CryoGridOutput from solution
out = @time solve(prob, Euler(), dt=120.0, saveat=24*3600.0, progress=true) |> CryoGridOutput;
# Plot it!
zs = [1,5,10,15,20:10:100...]
cg = Plots.cgrad(:copper,rev=true);
plot(ustrip.(out.H[Z(zs)]), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Enthalpy", leg=false, dpi=150)
plot(ustrip.(out.T[Z(zs)]), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)
