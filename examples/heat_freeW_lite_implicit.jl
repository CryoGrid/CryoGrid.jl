# # [Fast heat conduction with CryoGridLite](@id example7)
# This example is very similar to [Example 1](@ref) but uses the
# fast implicit CryoGridLite solver of Langer et al. 2023.

# Make sure to explicitly import the `LiteImplicit` submodule which has
# the relevant solver types.
using CryoGrid
using CryoGrid.LiteImplicit

# Load forcings and build stratigraphy like before.
forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_MkL3_CCSM4_long_term);
forcings = Base.rename(forcings, :Ptot => :precip)
z_top = -2.0u"m"
z_bot = 1000.0u"m"
upperbc = WaterHeatBC(
    SurfaceWaterBalance(forcings),
    TemperatureBC(forcings.Tair, NFactor()),
)
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
soil_layers = map(para -> Ground(para; heat, water), soilprofile)
strat = @Stratigraphy(
    z_top => Top(upperbc),
    z_top => Snowpack(Snow.LiteGridded(); heat),
    soil_layers...,
    z_bot => Bottom(GeothermalHeatFlux(0.053u"W/m^2"))
);
modelgrid = Grid(vcat(z_top:0.02u"m":-0.02u"m", CryoGrid.Presets.DefaultGrid_2cm))
tile = Tile(strat, modelgrid, ssinit);

# Since the solver can take daily timesteps, we can easily specify longer simulation time spans at minimal cost.
# Here we specify a time span of 30 years.
tspan = (DateTime(1990,10,1), DateTime(2020,10,1))
u0, du0 = initialcondition!(tile, tspan);
prob = CryoGridProblem(tile, u0, tspan, saveat=24*3600.0, savevars=(:T,:θw,:θwi,:kc,:ρsn))
sol = @time solve(prob, LiteImplicitEuler())
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
