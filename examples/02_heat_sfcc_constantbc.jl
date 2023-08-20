# # [Soil heat with SFCC and constant BCs](@id example2)
# In this example, we use the preset `SoilHeatTile` to construct
# a `Tile` consisting of a soil column with heat conduction
# and zero-flux boundary conditions. This is a useful test case for
# checking energy conservation since we can guarantee that no energy
# is being added or removed at the boundaries.

using CryoGrid

# Select default grid and initial temperature profile.
grid = CryoGrid.Presets.DefaultGrid_2cm
tempprofile = TemperatureProfile(
    0.0u"m" => -5.0u"°C",
    10.0u"m" => -1.0u"°C",
    100.0u"m" => 0.0u"°C",
);

# Here we use a simple single layer model with default soil parameters (50% porosity, no organic).
soilprofile = SoilProfile(
    0.0u"m" => MineralOrganic()
);

# Here we specify the soil freezing characteristic curve (SFCC) formulation of Painter and Karra (2014).
# The van Genuchten parameters `α=0.5` and `n=1.8` correspond to a silty soil.
sfcc = PainterKarra(swrc=VanGenuchten(α=0.5, n=1.8));

# Enthalpy form of the heat transfer operator (i.e. prognostic :H). In this case, this is equivalent to
# the shorthand `SoilHeatTile(:H, ...)`. However, it's worth demonstrating how the operator can be explicitly
# specified.
heatop = Heat.EnthalpyForm(SFCCNewtonSolver())
initT = initializer(:T, tempprofile)
tile = CryoGrid.Presets.SoilHeatTile(
    heatop,
    ConstantBC(HeatBalance, CryoGrid.Neumann, 0.0u"W/m^2"),
    ConstantBC(HeatBalance, CryoGrid.Neumann, 0.0u"W/m^2"),
    soilprofile,
    initT;
    grid, 
    freezecurve=sfcc,
    cachetype=CryoGrid.Numerics.ArrayCache,
);

# Define the simulation time span.
tspan = (DateTime(2010,1,1),DateTime(2010,12,31))
u0, du0 = initialcondition!(tile, tspan);

# Construct and solve the `CryoGridProblem`:
prob = CryoGridProblem(tile, u0, tspan, saveat=3600.0, savevars=(:T,))
@info "Running model"
sol = @time solve(prob, CGEuler(), dt=120.0, progress=true);
out = CryoGridOutput(sol)

# Compute total energy balance error.
Htot = Diagnostics.integrate(out.H, grid)
mass_balance_error = Htot[end] - Htot[1]

# Plot it!
import Plots
zs = [1,51,101]u"cm"
Diagnostics.plot_at_depths(:T, out, zs, ylabel="Temperature", leg=false, size=(800,500), dpi=150)

# Plot the energy balance error over time.
Plots.plot(uconvert.(u"MJ", Htot .- Htot[1]), title="Energy balance error")
