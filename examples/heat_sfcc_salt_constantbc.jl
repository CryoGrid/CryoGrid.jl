using CryoGrid
using Plots

grid = CryoGrid.Presets.DefaultGrid_2cm
tempprofile = TemperatureProfile(
    0.0u"m" => -1.0u"°C",
    20.0u"m" => -5.0u"°C",
    500.0u"m" => 0.0u"°C",
    1000.0u"m" => 10.2u"°C",
)
initT = initializer(:T, tempprofile)
initpor = SedimentCompactionInitializer(porosityZero=0.4)
sfcc = DallAmicoSalt(swrc=VanGenuchten(α=0.05, n=1.8))
sediment = MarineSediment(heat=HeatBalance(:T, freezecurve=sfcc), salt=SaltMassBalance())
strat = @Stratigraphy(
    0.0u"m" => Top(ConstantTemperature(0.0), SaltGradient(benthicSalt=900.0, surfaceState=0)),
    0.0u"m" => :sediment => sediment,
    1000.0u"m" => Bottom(GeothermalHeatFlux(0.053u"W/m^2"))
);
tile = Tile(strat, grid, initT, initpor)
# define time span
tspan = (DateTime(2010,1,1),DateTime(2010,1,2))
u0, du0 = initialcondition!(tile, tspan)
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(tile, u0, tspan, saveat=900.0, savevars=(:T,))
@info "Running model"
out = @time solve(prob, SSPRK43(), reltol=1e-8, saveat=900.0, progress=true) |> CryoGridOutput;
# Plot it!
zs = [5,10,15,20,25,30,40,50,100]u"cm"
cg = Plots.cgrad(:copper,rev=true);
plot(out.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)
# plot(out.θw[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)
Htot = Diagnostics.integrate(out.H, grid)
plot(uconvert.(u"MJ", Htot .- Htot[1]), title="Energy balance error")
# compute final energy balance error
mass_balance_error = Htot[end] - Htot[1]
