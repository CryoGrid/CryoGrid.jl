using CryoGrid
using NonlinearSolve

grid = CryoGrid.Presets.DefaultGrid_2cm
tempprofile = TemperatureProfile(
    0.0u"m" => -5.0u"°C",
    10.0u"m" => -1.0u"°C",
    100.0u"m" => 0.0u"°C",
)
soilprofile = SoilProfile(
    0.0u"m" => MineralOrganic()
)
initT = initializer(:T, tempprofile)
sfcc = PainterKarra(swrc=VanGenuchten(α=0.5, n=1.8))
heatop = Heat.EnthalpyForm(SFCCPreSolver())
tile = CryoGrid.Presets.SoilHeatTile(
    heatop,
    ConstantBC(HeatBalance, CryoGrid.Neumann, 0.0u"W/m^2"),
    ConstantBC(HeatBalance, CryoGrid.Neumann, 0.0u"W/m^2"),
    soilprofile,
    initT;
    grid, 
    freezecurve=sfcc
)
# define time span
tspan = (DateTime(2010,1,1),DateTime(2010,3,31))
u0, du0 = initialcondition!(tile, tspan)
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(tile, u0, tspan, saveat=900.0, savevars=(:T,), step_limiter=nothing)
@info "Running model"
out = @time solve(prob, Euler(), dt=120.0, saveat=900.0, progress=true) |> CryoGridOutput;

# check mass conservation
Htot = Diagnostics.integrate(out.H, grid)
# compute final energy balance error
mass_balance_error = Htot[end] - Htot[1]

# Plot it!
import Plots

zs = [1,3,5,7,11,31,51,101]u"cm"
Diagnostics.plot_at_depths(:T, out, zs, ylabel="Temperature", leg=false, size=(800,500), dpi=150)

# energy balance error over time
Plots.plot(uconvert.(u"MJ", Htot .- Htot[1]), title="Energy balance error")
