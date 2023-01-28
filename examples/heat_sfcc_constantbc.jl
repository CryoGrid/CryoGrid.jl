using CryoGrid
using Plots

grid = CryoGrid.Presets.DefaultGrid_2cm
tempprofile = TemperatureProfile(
    0.0u"m" => -1.0u"°C",
    100.0u"m" => 1.0u"°C",
)
soilprofile = SoilProfile(
    0.0u"m" => HomogeneousMixture()
)
initT = initializer(:T, tempprofile)
sfcc = SFCC(DallAmico(swrc=VanGenuchten(α=0.05, n=1.8)))
tile = CryoGrid.Presets.SoilHeatTile(
    :T,
    ConstantBC(HeatBalance, CryoGrid.Neumann, 0.0u"W/m^2"),
    ConstantBC(HeatBalance, CryoGrid.Neumann, 0.0u"W/m^2"),
    soilprofile,
    initT;
    grid=grid, 
    freezecurve=sfcc
)
# define time span
tspan = (DateTime(2010,1,1),DateTime(2010,12,31))
u0, du0 = initialcondition!(tile, tspan)
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(tile, u0, tspan, savevars=(:T,:H,:jH,:∂H∂t,:θw), step_limiter=nothing, saveat=900.0)
@info "Running model"
out = @time solve(prob, ImplicitEuler(), abstol=1e-8, reltol=1e-10, saveat=900.0, progress=true) |> CryoGridOutput;
# Plot it!
zs = [5,10,15,20,25,30,40,50,100]u"cm"
cg = Plots.cgrad(:copper,rev=true);
plot(out.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)
# plot(out.θw[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)
Htot = Diagnostics.integrate(out.H, grid)
plot(uconvert.(u"MJ", Htot .- Htot[1]), title="Energy balance error")
@show Htot[end] - Htot[1]
