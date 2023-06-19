# # Example 9
# ## Coupled heat and salt diffusion on salty soil column
# In this example, we construct a `Tile` from a Stratigraphy of three "salty soil" layers
# which include coupled salt/heat diffusion. Note that this is currently only supported for
# the temperature-based form of the heat equation.

using CryoGrid

initT = initializer(:T, -2.0)
initsalt = initializer(:c, 0.0)
initpor = SedimentCompactionInitializer(porosityZero=0.4)
sfcc = DallAmicoSalt(swrc=VanGenuchten(α=4.06, n=2.03))
upperbc = SaltHeatBC(SaltGradient(benthicSalt=900.0, surfaceState=0), ConstantTemperature(0.0))
strat = @Stratigraphy(
    0.0u"m" => Top(upperbc),
    0.0u"m" => :sediment => SaltySoil(heat=HeatBalance(:T, freezecurve=sfcc), salt=SaltMassBalance()),
    1000.0u"m" => Bottom(GeothermalHeatFlux(0.053u"W/m^2"))
);
modelgrid = CryoGrid.Presets.DefaultGrid_10cm
tile = Tile(strat, modelgrid, initT, initsalt, initpor)
# define time span
tspan = (DateTime(1990,1,1),DateTime(2000,12,31))
u0, du0 = initialcondition!(tile, tspan)
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(tile, u0, tspan, saveat=24*3600.0, savevars=(:T,:θw,:k,:dₛ))
@info "Running model"
integrator = init(prob, Euler(), dt=60.0, saveat=24*3600.0, progress=true);
# step forawrd 24 hours
@time step!(integrator, 24*3600)
# run to end of tspan
@time for i in integrator end
out = CryoGridOutput(integrator.sol)

# Plot it!
import Plots

zs = [5,15,25,35,45,55,65,75]u"cm"
cg = Plots.cgrad(:copper,rev=true);
total_salt = out.c.*out.θw
Plots.plot(total_salt[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Salt concentration", leg=false, size=(800,500), dpi=150)
Plots.plot(out.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)
Plots.plot(out.c[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Salt concentration", leg=false, size=(800,500), dpi=150)
# Plots.plot(out.θw[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Liquid water content", leg=false, size=(800,500), dpi=150)
