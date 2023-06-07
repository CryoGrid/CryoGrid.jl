using CryoGrid

grid = CryoGrid.Presets.DefaultGrid_2cm
tempprofile = TemperatureProfile(
    0.0u"m" => -1.0u"°C",
    100.0u"m" => 1.0u"°C",
)
soilprofile = SoilProfile(
    0.0u"m" => MineralOrganic()
)
initT = initializer(:T, tempprofile)
sfcc = DallAmico(swrc=VanGenuchten(α=0.05, n=1.8))
tile = CryoGrid.Presets.SoilHeatTile(
    :H,
    ConstantBC(HeatBalance, CryoGrid.Neumann, 0.0u"W/m^2"),
    ConstantBC(HeatBalance, CryoGrid.Neumann, 0.0u"W/m^2"),
    soilprofile,
    initT;
    grid, 
    freezecurve=sfcc
)
# define time span
tspan = (DateTime(2010,1,1),DateTime(2010,12,31))
u0, du0 = initialcondition!(tile, tspan)
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(tile, u0, tspan, saveat=900.0, savevars=(:T,))
@info "Running model"
out = @time solve(prob, SSPRK43(), reltol=1e-8, saveat=900.0, progress=true) |> CryoGridOutput;

# check mass conservation
Htot = Diagnostics.integrate(out.H, grid)
# compute final energy balance error
mass_balance_error = Htot[end] - Htot[1]

# Plot it!
import Plots

zs = [5,10,15,20,25,30,40,50,100]u"cm"
Diagnostics.plot_at_depths(:T, out, zs, ylabel="Temperature", leg=false, size=(800,500), dpi=150)

# energy balance error over time
Plots.plot(uconvert.(u"MJ", Htot .- Htot[1]), title="Energy balance error")
