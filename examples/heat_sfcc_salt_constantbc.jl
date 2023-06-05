using CryoGrid
using Plots

tempprofile = TemperatureProfile(
    0.0u"m" => -1.0u"°C",
    20.0u"m" => -5.0u"°C",
    500.0u"m" => 0.0u"°C",
    1000.0u"m" => 10.2u"°C",
)
initT = initializer(:T, tempprofile)
initsalt = initializer(:c, 800.0)
initpor = SedimentCompactionInitializer(porosityZero=0.4)
sfcc = DallAmicoSalt(swrc=VanGenuchten(α=4.06, n=2.03))
strat = @Stratigraphy(
    0.0u"m" => Top(ConstantTemperature(0.0), SaltGradient(benthicSalt=900.0, surfaceState=0)),
    0.0u"m" => :sediment => SaltySoil(heat=HeatBalance(:T, freezecurve=sfcc), salt=SaltMassBalance()),
    1000.0u"m" => Bottom(GeothermalHeatFlux(0.053u"W/m^2"))
);
modelgrid = CryoGrid.Presets.DefaultGrid_10cm
tile = Tile(strat, modelgrid, initT, initsalt, initpor)
# define time span
tspan = (DateTime(2000,1,1),DateTime(2000,12,31))
u0, du0 = initialcondition!(tile, tspan)
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(tile, u0, tspan, saveat=3600.0, savevars=(:T,:k,:dₛ))
@info "Running model"
integrator = init(prob, Euler(), dt=60.0, saveat=3600.0, progress=true);
# step forawrd 24 hours
@time step!(integrator, 24*3600)
# run to end of tspan
@time for i in integrator end
out = CryoGridOutput(integrator.sol)
# Plot it!
zs = [5,15,25,35,45,55,100,150,200]u"cm"
cg = Plots.cgrad(:copper,rev=true);
plot(out.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)
plot(out.c[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Salt concentration", leg=false, size=(800,500), dpi=150)
