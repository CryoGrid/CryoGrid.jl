# # [Soil heat with SFCC and constant BCs](@id example2)
# In this example, we use the preset `SoilHeatTile` to construct
# a `Tile` consisting of a soil column with heat conduction
# and zero-flux boundary conditions. This is a useful test case for
# checking energy conservation since we can guarantee that no energy
# is being added or removed at the boundaries.

using CryoGrid
using OrdinaryDiffEq

# Select default grid and initial temperature profile.
grid = CryoGrid.Presets.DefaultGrid_5cm
tempprofile = TemperatureProfile(
    0.0u"m" => -10.0u"°C",
    1000.0u"m" => 10.0u"°C",
);

# Here we use a simple single layer model with default soil parameters (50% porosity, no organic).
soilprofile = SoilProfile(
    0.0u"m" => MineralOrganic()
);

# Here we specify the soil freezing characteristic curve (SFCC) formulation of Painter and Karra (2014).
# The van Genuchten parameters `α=0.1` and `n=1.8` correspond to a silty soil.
sfcc = PainterKarra(swrc=VanGenuchten(α=0.1, n=1.8));

import Plots
Plots.plot(-2.0u"°C":0.01u"K":0.0u"°C", sfcc)

# Enthalpy form of the heat transfer operator (i.e. prognostic :H). In this case, this is equivalent to
# the shorthand `SoilHeatTile(:H, ...)`. However, it's worth demonstrating how the operator can be explicitly
# specified.
heatop = Heat.MOLEnthalpy(SFCCPreSolver())
initT = initializer(:T, tempprofile)
tile = CryoGrid.Presets.SoilHeatTile(
    heatop,
    ## 1 W/m^2 in and out, i.e. net zero flux
    ConstantBC(HeatBalance, CryoGrid.Neumann, 1.0u"W/m^2"),
    ConstantBC(HeatBalance, CryoGrid.Neumann, 1.0u"W/m^2"),
    soilprofile,
    initT;
    grid, 
    freezecurve=sfcc,
);

# Define the simulation time span.
tspan = (DateTime(2010,1,1),DateTime(2010,12,31))
u0, du0 = initialcondition!(tile, tspan);

# Construct and solve the `CryoGridProblem` using the 2nd order, implicit Trapezoid (Crank-Nicolson) solver from OrdinaryDiffEq
prob = CryoGridProblem(tile, u0, tspan, saveat=3600.0, savevars=(:T,:H), step_limiter=nothing)
sol = @time solve(prob, Trapezoid(), reltol=1e-8);
out = CryoGridOutput(sol)

# Compute total energy balance error.
Htot = Diagnostics.integrate(out.H, grid)
mass_balance_error = Htot[end] - Htot[1]

# Plot it!
import Plots
zs = [1,51,101]u"cm"
Diagnostics.plot_at_depths(:T, out, zs, ylabel="Temperature", leg=false, size=(800,500), dpi=150)

# Plot the energy balance error over time.
Plots.plot(uconvert.(u"MJ", Htot .- Htot[1]), title="", ylabel="Energy balance error")
