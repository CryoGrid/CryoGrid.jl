# # [Soil heat w/ SEB, snow cover, and bucket water scheme](@id example6)
# In this example, we construct a `Tile` consisting of a soil column with (i) heat conduction
# forced by the surface energy balance (SEB), (ii) a bulk snow scheme, and
# (iii) a bucket hydrology scheme.

# For this example, we need to use an OrdinaryDiffEq integrator.
using CryoGrid
using OrdinaryDiffEq

# First, load the forcings and construct the Tile.
modelgrid = CryoGrid.Presets.DefaultGrid_2cm;
soilprofile = SoilProfile(
    0.0u"m" => MineralOrganic(por=0.80,sat=0.8,org=0.75),
    0.1u"m" => MineralOrganic(por=0.80,sat=0.9,org=0.25),
    0.4u"m" => MineralOrganic(por=0.55,sat=1.0,org=0.25),
    3.0u"m" => MineralOrganic(por=0.50,sat=1.0,org=0.0),
    10.0u"m" => MineralOrganic(por=0.30,sat=1.0,org=0.0),
);
## mid-winter temperature profile
tempprofile = CryoGrid.Presets.SamoylovDefault.tempprofile
forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044);
tempprofile = CryoGrid.Presets.SamoylovDefault.tempprofile
initT = initializer(:T, tempprofile)
z = 2.0u"m"; # height [m] for which the forcing variables (Temp, humidity, wind, pressure) are provided
seb = SurfaceEnergyBalance(forcings, z)
swb = SurfaceWaterBalance(forcings)
upperbc = WaterHeatBC(swb, seb)
heat = HeatBalance(:H)
water = WaterBalance(BucketScheme(), DampedET())
## build stratigraphy
strat = @Stratigraphy(
    -z => Top(upperbc), 
    -z => Snowpack(heat=HeatBalance(), water=water),
    soilprofile[1].depth => Ground(soilprofile[1].value; heat, water),
    soilprofile[2].depth => Ground(soilprofile[2].value; heat, water),
    soilprofile[3].depth => Ground(soilprofile[3].value; heat, water),
    soilprofile[4].depth => Ground(soilprofile[4].value; heat, water),
    soilprofile[5].depth => Ground(soilprofile[5].value; heat, water),
    1000.0u"m" => Bottom(GeothermalHeatFlux(0.053u"J/s/m^2")),
);
## create Tile
tile = Tile(strat, modelgrid, initT);

# Set up the problem and solve it!
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
integrator = init(prob, Euler(), dt=60.0)
## step forwards 24 hours and check for NaN/Inf values
@time step!(integrator, 24*3600)
@assert all(isfinite.(integrator.u))
## iterate over remaining timespan at fixed points using `TimeChoiceIterator`
@time for (u,t) in TimeChoiceIterator(integrator, convert_t.(tspan[1]:Day(1):tspan[end]))
    @assert isfinite(getstate(:top, integrator).Qg[1])
    @info "Current t=$(Date(convert_t(t))), dt=$(integrator.dt)"
end
out = CryoGridOutput(integrator.sol)

# Plot it!
import Plots
zs = [1,5,10,15,20,25,30,40,50,100,150,200,500,1000]u"cm"
cg = Plots.cgrad(:copper,rev=true);
Plots.plot(ustrip.(out.T[Z(Near(zs))]), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)

# Saturation:
Plots.plot(ustrip.(out.sat[Z(Near(zs))]), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Soil saturation", leg=false, size=(800,500), dpi=150)

# Snow depth:
Plots.plot(ustrip.(out.snowpack.dsn), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Snow depth", leg=false, size=(800,500), dpi=150)

# Integrated ground heat flux:
Plots.plot(ustrip.(cumsum(out.top.Qg, dims=2)), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Integrated ground heat flux", leg=false, size=(800,500), dpi=150)
