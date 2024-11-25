# # [Computing parameter sensitivities with autodiff](@id example10)
# This example demonstrates how to parameterize and differentiate a simulation with
# two parameters (summer and winter n-factors) using forward-mode automatic simulation.
#
# TODO: add more detail/background

# Set up forcings and boundary conditions similarly to other examples:
using CryoGrid
forcings = loadforcings(CryoGrid.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044);
soilprofile, tempprofile = CryoGrid.Presets.SamoylovDefault
grid = CryoGrid.Presets.DefaultGrid_5cm
initT = initializer(:T, tempprofile)
tile = CryoGrid.Presets.SoilHeatTile(
    :T,
    TemperatureBC(Input(:Tair), NFactor(nf=Param(0.5), nt=Param(0.9))),
    GeothermalHeatFlux(0.053u"W/m^2"),
    soilprofile,
    forcings,
    initT;
    grid=grid
)
tspan = (DateTime(2010,10,1),DateTime(2010,10,2))
u0, du0 = @time initialcondition!(tile, tspan);

# Collect model parameters
p = CryoGrid.parameters(tile)

# Create the `CryoGridProblem`.
prob = CryoGridProblem(tile, u0, tspan, p, saveat=3600.0);

# Solve the forward problem with default parameter settings:
sol = @time solve(prob)
out = CryoGridOutput(sol)

# ForwardDiff provides tools for forward-mode automatic differentiation.
using ForwardDiff
using SciMLSensitivity
using Zygote

# Define a "loss" function; here we'll just take the mean over the final temperature field.
using Statistics
function loss(prob::CryoGridProblem, p)
    # local u0, _ = initialcondition!(tile, tspan, p)
    # local prob = CryoGridProblem(tile, u0, tspan, p, saveat=24*3600.0)
    newprob = remake(prob, p=p)
    sensealg = InterpolatingAdjoint(autojacvec=true, checkpointing=true)
    newsol = solve(newprob, Euler(), dt=300.0, sensealg=sensealg);
    newout = CryoGridOutput(newsol)
    return mean(ustrip.(newout.T[:,end]))
end

# Compute gradient with forward diff:
pvec = vec(p)
grad = @time ForwardDiff.gradient(pᵢ -> loss(prob, pᵢ), pvec)
grad = @time Zygote.gradient(pᵢ -> loss(prob, pᵢ), pvec)
