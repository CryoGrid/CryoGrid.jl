# # [Computing parameter sensitivities with autodiff](@id example10)
# This example demonstrates how to parameterize and differentiate a simulation with
# two parameters (summer and winter n-factors) using forward-mode automatic simulation.
#
# TODO: add more detail/background

# Set up forcings and boundary conditions similarly to other examples:
using CryoGrid
forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044);
soilprofile, tempprofile = CryoGrid.Presets.SamoylovDefault
grid = CryoGrid.Presets.DefaultGrid_5cm
initT = initializer(:T, tempprofile)
tile = CryoGrid.Presets.SoilHeatTile(
    :T,
    TemperatureGradient(forcings.Tair, NFactor(nf=Param(0.5), nt=Param(0.9))),
    GeothermalHeatFlux(0.053u"W/m^2"),
    soilprofile,
    initT;
    freezecurve=PainterKarra(),
    grid=grid
)
tspan = (DateTime(2010,10,30),DateTime(2011,10,30))
u0, du0 = initialcondition!(tile, tspan);

# Collect model parameters
p = CryoGrid.parameters(tile)

# Create the `CryoGridProblem`.
prob = CryoGridProblem(tile, u0, tspan, p, saveat=24*3600.0);

# Define a "loss" function; here we'll just take the mean over the final temperature field.
using Statistics
function loss(p)
    local u0, _ = initialcondition!(tile, tspan, p)
    local prob = CryoGridProblem(tile, u0, tspan, p, saveat=24*3600.0)
    local sol = solve(prob, CGEuler());
    local out = CryoGridOutput(sol)
    return mean(ustrip.(out.T[:,end]))
end

# ForwardDiff provides tools for forward-mode automatic differentiation.
using ForwardDiff
pvec = vec(p)

# Compute gradient with forward diff:
grad = ForwardDiff.gradient(loss, pvec)
