# # [Computing parameter sensitivities with autodiff](@id example10)
# This example demonstrates how to parameterize and differentiate a simulation with
# two parameters (summer and winter n-factors) using forward-mode automatic simulation.
#
# TODO: add more detail/background
using CryoGrid

# Set up forcings and boundary conditions similarly to other examples:
forcings = loadforcings(CryoGrid.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044);
soilprofile, tempprofile = CryoGrid.SamoylovDefault
freezecurve = PainterKarra(swrc=VanGenuchten())
soilprofile = SoilProfile(0.0u"m" => SimpleSoil(; freezecurve))
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
tspan = (DateTime(2010,9,1),DateTime(2011,10,1))
u0, du0 = @time initialcondition!(tile, tspan);

# We can retrieve the parameters of the system from `tile`:
para = CryoGrid.parameters(tile)

# Create the `CryoGridProblem`.
prob = CryoGridProblem(tile, u0, tspan, saveat=3600.0)

# Solve the forward problem with default parameter settings:
sol = @time solve(prob)
out = CryoGridOutput(sol)

# Import relevant packages for automatic differentiation.
using ForwardDiff
using SciMLSensitivity
using Zygote

# Define a "loss" function; here we'll just take the mean over the final temperature field.
using OrdinaryDiffEq
using Statistics
function loss(prob::CryoGridProblem, p)
    newprob = remake(prob, p=p)
    ## Here we specify the sensitivity algorithm. Note that this is only
    ## necessary for reverse-mode autodiff with Zygote.
    ## autojacvec = true uses ForwardDiff to calculate the jacobian;
    ## enabling checkpointing (theoretically) reduces the memory cost of the backwards pass.
    sensealg = InterpolatingAdjoint(autojacvec=true, checkpointing=true)
    newsol = solve(newprob, Euler(), dt=300.0, sensealg=sensealg);
    newout = CryoGridOutput(newsol)
    return mean(ustrip.(newout.T[:,end]))
end

# Compute gradient with forward-mode autodiff:
pvec = vec(prob.p)
fd_grad = @time ForwardDiff.gradient(pᵢ -> loss(prob, pᵢ), pvec)

## We can also try with reverse-mode autodiff. This is generally slower for smaller numbers
## of parmaeters (<100) but could be worthwhile for model configurations with high-dimensional
## parameterizations.
## zy_grad = @time Zygote.gradient(pᵢ -> loss(prob, pᵢ), pvec)
## @assert maximum(abs.(fd_grad .- zy_grad)) .< 1e-6 "Forward and reverse gradients don't match!"
