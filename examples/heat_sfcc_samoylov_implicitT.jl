# This example demonstrates how to simulate two-phase heat conduction without water flow
# and using the temperature-based formulation of the heat equation.
using CryoGrid
using OrdinaryDiffEq

import Plots

# Here we define the boundary conditions and stratigraphy.
forcings = loadforcings(CryoGrid.Forcings.Samoylov_ERA_MkL3_CCSM4_long_term);
upperbc = TemperatureBC(Input(:Tair), NFactor(0.65,0.9))
ssinit = ThermalSteadyStateInit(T0=-15.0u"°C");

# Create and plot the SFCC for sandy soil that we will use.
sfcc = PainterKarra(swrc=VanGenuchten("sand"))
Plots.plot(-0.5u"°C":0.001u"K":0.0u"°C", sfcc)

# Define soil stratigraphy.
soilprofile = SoilProfile(
    0.0u"m" => SimpleSoil(; por=0.80, org=0.75, freezecurve=sfcc),
    0.1u"m" => SimpleSoil(; por=0.80, org=0.25, freezecurve=sfcc),
    0.4u"m" => SimpleSoil(; por=0.55, org=0.25, freezecurve=sfcc),
    3.0u"m" => SimpleSoil(; por=0.50, org=0.0, freezecurve=sfcc),
    10.0u"m" => SimpleSoil(; por=0.30, org=0.0, freezecurve=sfcc),
)

# Set up model using the temperature formulation of the heat equation.
# This is also sometimes referred to as the "apparent heat capacity method".
heat = HeatBalance(:T)
soil_layers = map(para -> Ground(para; heat), soilprofile);
modelgrid = CryoGrid.DefaultGrid_2cm
strat = Stratigraphy(
    0.0u"m" => Top(upperbc),
    soil_layers,
    modelgrid[end] => Bottom(GeothermalHeatFlux(0.053u"W/m^2"))
);
tile = Tile(strat, modelgrid, forcings, ssinit);

tspan = (DateTime(2000,10,1), DateTime(2002,10,1))
u0, du0 = initialcondition!(tile, tspan);
prob = CryoGridProblem(tile, u0, tspan, saveat=24*3600, savevars=(:T,:θw,:∂H∂T), step_limiter=nothing)
integrator = init(prob, ImplicitEuler(autodiff=false))
# test one step
@time step!(integrator)
last_t = integrator.t
for i in integrator
    if integrator.t - last_t > 24*3600
        println("t=$(convert_t(integrator.t)), dt=$(integrator.dt)")
        last_t = integrator.t
    end
end
out = CryoGridOutput(integrator.sol)

# Plot the results!
zs = [1,5,10,15,20,25,30,40,50,100]u"cm"
cg = Plots.cgrad(:copper,rev=true);
Plots.plot(out.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", title="", leg=false, dpi=150)
