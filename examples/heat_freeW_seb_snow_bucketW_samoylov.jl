using CryoGrid
using Dates
using Plots

# Custom grid;
# Similar to Presets.DefaultGrid_2cm but includes one large cell above the surface for snow.
const gridvals = vcat([
    -2.0,
    0:0.02:2...,
    2.05:0.05:4.0...,
	4.1:0.1:10...,
    10.2:0.2:20...,
    21:1:30...,
	35:5:50...,
    60:10:100...,
    200:100:1000...]...
)u"m"
modelgrid = Grid(gridvals);
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
    :pressure => u"Pa",
    :wind => u"m/s",
    :Lin => u"W/m^2",
    :Sin => u"W/m^2",
    :snowfall => u"mm/d",
    :rainfall => u"mm/d",
);
Tair = TimeSeriesForcing(forcings.data.Tair, forcings.timestamps, :Tair);
pr   = TimeSeriesForcing(forcings.data.pressure, forcings.timestamps, :pr);
q    = TimeSeriesForcing(forcings.data.q, forcings.timestamps, :q);
wind = TimeSeriesForcing(forcings.data.wind, forcings.timestamps, :wind);
Lin  = TimeSeriesForcing(forcings.data.Lin, forcings.timestamps, :Lin);
Sin  = TimeSeriesForcing(forcings.data.Sin, forcings.timestamps, :Sin);
snowfall = TimeSeriesForcing(uconvert.(u"m/s", forcings.data.snowfall).*1.0, forcings.timestamps, :sf)
rainfall = TimeSeriesForcing(uconvert.(u"m/s", forcings.data.rainfall).*1.0, forcings.timestamps, :rf)
z = 2.;    # height [m] for which the forcing variables (Temp, humidity, wind, pressure) are provided
soilprofile, tempprofile = CryoGrid.Presets.SamoylovDefault
initT = initializer(:T, tempprofile)
# initialize saturation to match soil profile
initsat = initializer(:sat, (l,p,state) -> state.sat .= l.para.sat)
seb = SurfaceEnergyBalance(Tair, pr, q, wind, Lin, Sin, z)
subsurface_layers = map(enumerate(soilprofile)) do (i, soil_i)
    name = Symbol(:soil, i)
    heat = HeatBalance(:H, freezecurve=PainterKarra())
    water = WaterBalance(BucketScheme(), DampedET())
    soil_i.depth => name => Soil(soil_i.value; heat, water)
end
# build stratigraphy
strat = @Stratigraphy(
    -z*u"m" => Top(seb, Rainfall(rainfall), Snowfall(snowfall)),
    0.0u"m" => :snowpack => Snowpack(heat=HeatBalance()),
    subsurface_layers...,
    1000.0u"m" => Bottom(GeothermalHeatFlux(0.053u"J/s/m^2")),
);
# create Tile
tile = Tile(strat, modelgrid, initT, initsat);

# define time span
tspan = (DateTime(2010,10,30), DateTime(2011,10,30))
# generate initial condition and set up CryoGridProblem
u0, du0 = initialcondition!(tile, tspan)
prob = CryoGridProblem(tile, u0, tspan, savevars=(:T,:top => (:Qh,:Qe,:Qg,),:snowpack => (:dsn,)), step_limiter=nothing, saveat=3*3600.0)
# initialize integrator
integrator = init(prob, Euler(), dt=120.0, progress=true)
# take one 24 hour step to check that everything is working
@time step!(integrator, 24*3600)
@assert all(isfinite.(integrator.u))
# iterate over remaining timespan
@time for (u,t) in TimeChoiceIterator(integrator, convert_t.(tspan[1]:Day(1):tspan[end]))
    # print once per day to track progress
    @show convert_t(t)
end
# build output from solution
out = CryoGridOutput(integrator.sol)
# Plot it!
zs = [1,5,10,15,20,25,30,40,50]u"cm"
cg = Plots.cgrad(:copper,rev=true);
plot(ustrip.(out.H[Z(Near(zs))]), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Enthalpy", leg=false, dpi=150)
plot(ustrip.(out.T[Z(Near(zs))]), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)
plot(ustrip.(out.sat[Z(Near(zs))]), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Soil saturation", leg=false, size=(800,500), dpi=150)
plot(ustrip.(out.snowpack.dsn), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Snow depth", leg=false, size=(800,500), dpi=150)
plot(ustrip.(cumsum(out.top.Qg, dims=2)), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Integrated ground heat flux", leg=false, size=(800,500), dpi=150)
