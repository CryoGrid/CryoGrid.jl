# # [Computing parameter sensitivities with autodiff](@id example10)
# This example demonstrates how to parameterize and differentiate a simulation with
# two parameters (summer and winter n-factors) using forward-mode automatic simulation.
#
# TODO: add more detail/background
using CryoGrid

# Set up forcings and boundary conditions similarly to other examples:
forcings = loadforcings(CryoGrid.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044);
soilprofile, tempprofile = CryoGrid.SamoylovDefault
grid = CryoGrid.DefaultGrid_5cm
initT = initializer(:T, tempprofile)
tile = CryoGrid.SoilHeatTile(
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

# Create the `CryoGridProblem`.
prob = CryoGridProblem(tile, u0, tspan, saveat=3600.0);

# Solve the forward problem with default parameter settings:
sol = @time solve(prob)
out = CryoGridOutput(sol)

# Import relevant packages for automatic differentiation.
using ForwardDiff
using SciMLSensitivity
using Zygote

# Define a "loss" function; here we'll just take the mean over the final temperature field.
using Statistics
function loss(prob::CryoGridProblem, p)
    newprob = remake(prob, p=p)
    # autojacvec = true uses ForwardDiff to calculate the jacobian;
    # enabling checkpointing (theroetically) reduces the memory cost of the backwards pass.
    sensealg = InterpolatingAdjoint(autojacvec=true, checkpointing=true)
    newsol = solve(newprob, Euler(), dt=300.0, sensealg=sensealg);
    newout = CryoGridOutput(newsol)
    return mean(ustrip.(newout.T[:,end]))
end

# Compute gradient with forward diff:
pvec = prob.p
fd_grad = @time ForwardDiff.gradient(pᵢ -> loss(prob, pᵢ), pvec)
zy_grad = @time Zygote.gradient(pᵢ -> loss(prob, pᵢ), pvec)
@assert maximum(abs.(fd_grad .- zy_grad)) .< 1e-4 "Forward and reverse gradients don't match!"
@show fd_grad
