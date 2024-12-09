# # [Fast heat conduction with CryoGridLite](@id example7)
# This example is very similar to [Example 1](@ref) but uses the
# fast implicit CryoGridLite solver of Langer et al. 2023.

# Make sure to explicitly import the `LiteImplicit` submodule which has
# the relevant solver types.
using CryoGrid
using CryoGrid.LiteImplicit

# First we prepare the forcing data. This data includes only air temperature and
# *total* precipitation, so we need to separate this into rainfall and snowfall
# using air temperature.
raw_forcings = loadforcings(CryoGrid.Forcings.Samoylov_ERA_MkL3_CCSM4_long_term);
Tair = raw_forcings.data.Tair
Ptot = uconvert.(u"m/s", raw_forcings.data.Ptot)
rainfall = Ptot.*(Tair .> 0u"°C")
snowfall = Ptot.*(Tair .<= 0u"°C")
forcings = rebuild(raw_forcings; Tair, rainfall, snowfall);

# Next we define the boundary conditions and stratigraphy.
z_top = -2.0u"m"
z_bot = 1000.0u"m"
upperbc = WaterHeatBC(
    SurfaceWaterBalance(Input(:rainfall), Input(:snowfall)),
    TemperatureBC(Input(:Tair), NFactor(0.5,0.9))
);
ssinit = ThermalSteadyStateInit(T0=-15.0u"°C")
heatop = Heat.EnthalpyImplicit()
soilprofile = SoilProfile(
    0.0u"m" => SimpleSoil(; por=0.80, org=0.75, freezecurve=FreeWater()),
    0.1u"m" => SimpleSoil(; por=0.80, org=0.25, freezecurve=FreeWater()),
    0.4u"m" => SimpleSoil(; por=0.55, org=0.25, freezecurve=FreeWater()),
    3.0u"m" => SimpleSoil(; por=0.50, org=0.0, freezecurve=FreeWater()),
    10.0u"m" => SimpleSoil(; por=0.30, org=0.0, freezecurve=FreeWater()),
)
heat = HeatBalance(heatop)
water = WaterBalance()
soil_layers = map(para -> Ground(para; heat, water), soilprofile);
strat = @Stratigraphy(
    z_top => Top(upperbc),
    z_top => Snowpack(Snow.LiteGridded(); heat),
    soil_layers...,
    z_bot => Bottom(GeothermalHeatFlux(0.053u"W/m^2"))
);
modelgrid = Grid(vcat(z_top:0.02u"m":-0.02u"m", CryoGrid.DefaultGrid_2cm))
tile = Tile(strat, modelgrid, forcings, ssinit);

# Since the solver can take daily timesteps, we can easily specify longer simulation time spans at minimal cost.
# Here we specify a time span of 20 years.
tspan = (DateTime(2000,10,1), DateTime(2020,10,1))
u0, du0 = initialcondition!(tile, tspan);
prob = CryoGridProblem(tile, u0, tspan, saveat=24*3600, savevars=(:T,:θw))
sol = @time solve(prob, LiteImplicitEuler(), dt=24*3600)
out = CryoGridOutput(sol)

# Plot the results!
import Plots
zs = [1,5,10,15,20,25,30,40,50,100]u"cm"
cg = Plots.cgrad(:copper,rev=true);
Plots.plot(out.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", title="", leg=false, dpi=150)
Plots.plot(out.snowpack.dsn)

# CryoGridLite can also be embedded into integrators from OrdinaryDiffEq.jl via the `NLCGLite` nonlinear solver interface.
# Note that these sovers generally will not be faster (in execution time) but may be more stable in some cases. Adaptive timestepping can be employed by
# removing the `adaptive=false` argument.
# using OrdinaryDiffEq
# sol2 = @time solve(prob, ImplicitEuler(nlsolve=NLCGLite()), adaptive=false, dt=24*3600.0, saveat=24*3600.0);
