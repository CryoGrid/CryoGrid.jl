# # [Soil heat with SFCC and constant BCs](@id example2)
# In this example, we use the preset `SoilHeatTile` to construct
# a `Tile` consisting of a soil column with heat conduction
# and zero-flux boundary conditions. This is a useful test case for
# checking energy conservation since we can guarantee that no energy
# is being added or removed at the boundaries.

using CryoGrid

# Select default grid with 2 cm near-surface spacing.
grid = CryoGrid.Presets.DefaultGrid_2cm
# Specify initial temperature and soil profiles.
tempprofile = TemperatureProfile(
    0.0u"m" => -5.0u"°C",
    10.0u"m" => -1.0u"°C",
    100.0u"m" => 0.0u"°C",
)
soilprofile = SoilProfile(
    0.0u"m" => MineralOrganic()
)
initT = initializer(:T, tempprofile)
# Here we specify the soil freezing characteristic curve (SFCC) formulation of Painter and Karra (2014).
# The van Genuchten parameters `α=0.5` and `n=1.8` correspond to a silty soil.
sfcc = PainterKarra(swrc=VanGenuchten(α=0.5, n=1.8))
# Enthalpy form of the heat transfer operator (i.e. prognostic :H). In this case, this is equivalent to
# the shorthand `SoilHeatTile(:H, ...)`. However, it's worth demonstrating how the operator can be explicitly
# specified.
heatop = Heat.EnthalpyForm(SFCCNewtonSolver())
tile = CryoGrid.Presets.SoilHeatTile(
    heatop,
    ConstantBC(HeatBalance, CryoGrid.Neumann, 0.0u"W/m^2"),
    ConstantBC(HeatBalance, CryoGrid.Neumann, 0.0u"W/m^2"),
    soilprofile,
    initT;
    grid, 
    freezecurve=sfcc,
    cachetype=CryoGrid.Numerics.ArrayCache,
)
# Define the simulation time span.
tspan = (DateTime(2010,1,1),DateTime(2010,12,31))
u0, du0 = initialcondition!(tile, tspan)
# Construct and solve the `CryoGridProblem`:
prob = CryoGridProblem(tile, u0, tspan, saveat=3600.0, savevars=(:T,))
@info "Running model"
out = @time solve(prob, Euler(), dt=120.0, saveat=3600.0, progress=true) |> CryoGridOutput;

# check mass conservation
Htot = Diagnostics.integrate(out.H, grid)
# compute final energy balance error
mass_balance_error = Htot[end] - Htot[1]

# Plot it!
import Plots

zs = [1,51,101]u"cm"
Diagnostics.plot_at_depths(:T, out, zs, ylabel="Temperature", leg=false, size=(800,500), dpi=150)

# energy balance error over time
Plots.plot(uconvert.(u"MJ", Htot .- Htot[1]), title="Energy balance error")
