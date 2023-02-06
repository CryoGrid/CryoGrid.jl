using CryoGrid
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
grid = CryoGrid.Presets.DefaultGrid_5cm
_, tempprofile = CryoGrid.Presets.SamoylovDefault
# soil profile: depth => (excess ice, natural porosity, saturation, organic fraction)
soilprofile = SoilProfile(
    0.0u"m" => HomogeneousMixture(por=0.80,sat=0.5,org=0.75), #(θwi=0.80,θm=0.05,θo=0.15,ϕ=0.80),
    0.1u"m" => HomogeneousMixture(por=0.80,sat=0.7,org=0.25), #(θwi=0.80,θm=0.15,θo=0.05,ϕ=0.80),
    0.4u"m" => HomogeneousMixture(por=0.55,sat=0.9,org=0.25), #(θwi=0.80,θm=0.15,θo=0.05,ϕ=0.55),
    3.0u"m" => HomogeneousMixture(por=0.50,sat=1.0,org=0.0), #(θwi=0.50,θm=0.50,θo=0.0,ϕ=0.50),
    10.0u"m" => HomogeneousMixture(por=0.30,sat=1.0,org=0.0), #(θwi=0.30,θm=0.70,θo=0.0,ϕ=0.30),
    100.0u"m" => HomogeneousMixture(por=0.10,sat=1.0,org=0.0),
);
initT = initializer(:T, tempprofile)
initsat = initializer(:sat, (l,p,state) -> state.sat .= l.para.sat)
# @Stratigraphy macro lets us list multiple subsurface layers
strat = @Stratigraphy(
    -2.0*u"m" => Top(TemperatureGradient(tair), Rainfall(pr)),
    soilprofile[1].depth => :soil1 => Soil(Coupled(WaterBalance(BucketScheme()), HeatBalance()), para=soilprofile[1].value),
    soilprofile[2].depth => :soil2 => Soil(Coupled(WaterBalance(BucketScheme()), HeatBalance()), para=soilprofile[2].value),
    soilprofile[3].depth => :soil3 => Soil(Coupled(WaterBalance(BucketScheme()), HeatBalance()), para=soilprofile[3].value),
    soilprofile[4].depth => :soil4 => Soil(Coupled(WaterBalance(BucketScheme()), HeatBalance()), para=soilprofile[4].value),
    soilprofile[5].depth => :soil5 => Soil(Coupled(WaterBalance(BucketScheme()), HeatBalance()), para=soilprofile[5].value),
    1000.0u"m" => Bottom(GeothermalHeatFlux(0.053u"J/s/m^2")),
);
tile = Tile(strat, grid, initT, initsat);
# define time span
tspan = (DateTime(2011,10,30),DateTime(2012,10,30))
u0, du0 = initialcondition!(tile, tspan)
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(tile, u0, tspan, savevars=(:T,:θw,:θwi), saveat=3*3600.0)
out = @time solve(prob, Euler(), dt=300.0, saveat=3*3600.0, progress=true) |> CryoGridOutput;
# check mass conservation
water_added = values(sum(pr[tspan[1]:Hour(3):tspan[2]].tarray.*(3*3600.0u"s")))[1]
water_mass = Diagnostics.integrate(out.θwi, tile.grid)
Δwater = water_mass[end] - water_mass[1]
# Plot it!
zs = [1,5,10,15,20,30,40,50,100,150,200]u"cm"
cg = Plots.cgrad(:copper,rev=true);
plot(out.H[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Enthalpy", leg=false, dpi=150)
plot(out.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)
plot(out.θw[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Unfrozen water content", leg=false, size=(800,500), dpi=150)
plot(out.sat[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Saturation", leg=false, size=(800,500), dpi=150)
