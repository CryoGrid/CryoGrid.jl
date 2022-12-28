using CryoGrid
using Plots

forcings = loadforcings(
    CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044,
    :Tair => u"°C",
    :rainfall => u"mm",
);
# use air temperature as upper boundary forcing;
tair = TimeSeriesForcing(forcings.data.Tair, forcings.timestamps, :Tair);
pr = TimeSeriesForcing(uconvert.(u"m/s", forcings.data.rainfall./3u"hr"), forcings.timestamps, :rainfall)
# "simple" heat conduction model w/ 5 cm grid spacing
grid = CryoGrid.Presets.DefaultGrid_5cm
_, tempprofile = CryoGrid.Presets.SamoylovDefault
# soil profile: depth => (excess ice, natural porosity, saturation, organic fraction)
soilprofile = SoilProfile(
    0.0u"m" => HomogeneousMixture(por=0.80,sat=0.5,org=0.75), #(θwi=0.80,θm=0.05,θo=0.15,ϕ=0.80),
    0.1u"m" => HomogeneousMixture(por=0.80,sat=0.7,org=0.25), #(θwi=0.80,θm=0.15,θo=0.05,ϕ=0.80),
    0.4u"m" => HomogeneousMixture(por=0.55,sat=0.9,org=0.25), #(θwi=0.80,θm=0.15,θo=0.05,ϕ=0.55),
    3.0u"m" => HomogeneousMixture(por=0.50,sat=1.0,org=0.0), #(θwi=0.50,θm=0.50,θo=0.0,ϕ=0.50),
    10.0u"m" => HomogeneousMixture(por=0.30,sat=1.0,org=0.0), #(θwi=0.30,θm=0.70,θo=0.0,ϕ=0.30),
    100.0u"m" => HomogeneousMixture(por=0.10,sat=0.1,org=0.0),
);
initT = initializer(:T, tempprofile)
initsat = initializer(:sat, (l,p,state) -> state.sat .= l.para.sat)
# @Stratigraphy macro lets us list multiple subsurface layers
strat = @Stratigraphy(
    -2.0*u"m" => Top(TemperatureGradient(tair), Rainfall(pr)),
    soilprofile[1].depth => :soil1 => Soil(Coupled(WaterBalance(BucketScheme()), HeatBalance(freezecurve=SFCC(PainterKarra()))), para=soilprofile[1].value),
    soilprofile[2].depth => :soil2 => Soil(Coupled(WaterBalance(BucketScheme()), HeatBalance(freezecurve=SFCC(PainterKarra()))), para=soilprofile[2].value),
    soilprofile[3].depth => :soil3 => Soil(Coupled(WaterBalance(BucketScheme()), HeatBalance(freezecurve=SFCC(PainterKarra()))), para=soilprofile[3].value),
    soilprofile[4].depth => :soil4 => Soil(Coupled(WaterBalance(BucketScheme()), HeatBalance(freezecurve=SFCC(PainterKarra()))), para=soilprofile[4].value),
    soilprofile[5].depth => :soil5 => Soil(Coupled(WaterBalance(BucketScheme()), HeatBalance(freezecurve=SFCC(PainterKarra()))), para=soilprofile[5].value),
    1000.0u"m" => Bottom(GeothermalHeatFlux(0.053u"J/s/m^2")),
);
grid = CryoGrid.Presets.DefaultGrid_2cm
tile = Tile(strat, grid, initT, initsat);
# define time span
tspan = (DateTime(2010,10,30),DateTime(2011,10,30))
u0, du0 = initialcondition!(tile, tspan)
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(tile, u0, tspan, savevars=(:T,))
out = @time solve(prob, Euler(), dt=900.0, saveat=24*3600.0, progress=true) |> CryoGridOutput;
# Plot it!
zs = [1:10...,20:10:100...]
cg = Plots.cgrad(:copper,rev=true);
plot(out.H[Z(zs)], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Enthalpy", leg=false, dpi=150)
plot(out.T[Z(zs)], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)
plot(out.sat[Z(1:2:20)], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Saturation", leg=false, size=(800,500), dpi=150)
