# # [Soil heat w/ SEB, snow cover, and bucket water scheme](@id example6)
# In this example, we construct a `Tile` consisting of a soil column with (i) heat conduction
# forced by the surface energy balance (SEB), (ii) a bulk snow scheme, and
# (iii) a bucket hydrology scheme.

using CryoGrid

modelgrid = CryoGrid.Presets.DefaultGrid_2cm;
## soil profile: depth => (excess ice, natural porosity, saturation, organic fraction)
soilprofile = SoilProfile(
    0.0u"m" => MineralOrganic(por=0.80,sat=1.0,org=0.75), #(θwi=0.80,θm=0.05,θo=0.15,ϕ=0.80),
    0.1u"m" => MineralOrganic(por=0.80,sat=1.0,org=0.25), #(θwi=0.80,θm=0.15,θo=0.05,ϕ=0.80),
    0.4u"m" => MineralOrganic(por=0.55,sat=1.0,org=0.25), #(θwi=0.80,θm=0.15,θo=0.05,ϕ=0.55),
    3.0u"m" => MineralOrganic(por=0.50,sat=1.0,org=0.0), #(θwi=0.50,θm=0.50,θo=0.0,ϕ=0.50),
    10.0u"m" => MineralOrganic(por=0.30,sat=1.0,org=0.0), #(θwi=0.30,θm=0.70,θo=0.0,ϕ=0.30),
);
## mid-winter temperature profile
tempprofile = CryoGrid.Presets.SamoylovDefault.tempprofile
forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044);
soilprofile, tempprofile = CryoGrid.Presets.SamoylovDefault
initT = initializer(:T, tempprofile)
## initialize saturation to match soil profile
initsat = initializer(:sat, (l,state) -> state.sat .= l.para.sat)
z = 2.0u"m"; # height [m] for which the forcing variables (Temp, humidity, wind, pressure) are provided
seb = SurfaceEnergyBalance(forcings.Tair, forcings.pressure, forcings.q, forcings.wind, forcings.Lin, forcings.Sin, z)
swb = SurfaceWaterBalance(rainfall=forcings.rainfall, snowfall=forcings.snowfall)
upperbc = WaterHeatBC(swb, seb)
heat = HeatBalance(:H, freezecurve=PainterKarra())
water = WaterBalance(BucketScheme(), DampedET())
## build stratigraphy
strat = @Stratigraphy(
    -z => Top(upperbc),
    -z => :snowpack => Snowpack(heat=HeatBalance()),
    soilprofile[1].depth => :soil1 => SimpleSoil(soilprofile[1].value; heat, water),
    soilprofile[2].depth => :soil2 => SimpleSoil(soilprofile[2].value; heat, water),
    soilprofile[3].depth => :soil3 => SimpleSoil(soilprofile[3].value; heat, water),
    soilprofile[4].depth => :soil4 => SimpleSoil(soilprofile[4].value; heat, water),
    soilprofile[5].depth => :soil5 => SimpleSoil(soilprofile[5].value; heat, water),
    1000.0u"m" => Bottom(GeothermalHeatFlux(0.053u"J/s/m^2")),
);
## create Tile
tile = Tile(strat, modelgrid, initT, initsat);

## define time span
tspan = (DateTime(2010,10,30), DateTime(2011,10,30))
## generate initial condition and set up CryoGridProblem
u0, du0 = initialcondition!(tile, tspan)
prob = CryoGridProblem(
    tile,
    u0,
    tspan,
    savevars=(:T,:jH,:top => (:Qh,:Qe,:Qg,),:snowpack => (:dsn,)),
    saveat=3*3600.0
)
## initialize integrator
integrator = init(prob, CGEuler(), dt=60.0)
## step forwards 24 hours and check for NaN/Inf values
@time step!(integrator, 24*3600)
@assert all(isfinite.(integrator.u))
## iterate over remaining timespan
@time for (u,t) in TimeChoiceIterator(integrator, convert_t.(tspan[1]:Day(1):tspan[end]))
    @assert isfinite(getstate(:top, integrator).Qg[1])
    @show Date(convert_t(t))
end
## build output from solution
out = CryoGridOutput(integrator.sol)

# Plot it!
import Plots

zs = [1,5,10,15,20,25,30,40,50]u"cm"
cg = Plots.cgrad(:copper,rev=true);
Plots.plot(ustrip.(out.T[Z(Near(zs))]), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)
Plots.plot(ustrip.(out.sat[Z(Near(zs))]), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Soil saturation", leg=false, size=(800,500), dpi=150)
Plots.plot(ustrip.(out.snowpack.dsn), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Snow depth", leg=false, size=(800,500), dpi=150)
Plots.plot(ustrip.(cumsum(out.top.Qg, dims=2)), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Integrated ground heat flux", leg=false, size=(800,500), dpi=150)
