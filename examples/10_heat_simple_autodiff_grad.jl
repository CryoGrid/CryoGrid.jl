# # [Computing parameter sensitivities with autodiff](@id example10)
# This example demonstrates how to parameterize and differentiate a simulation with
# two parameters (summer and winter n-factors) using forward-mode automatic simulation.
#
# TODO: add more detail/background

using CryoGrid
using Statistics

forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044, :Tair => u"°C");
soilprofile, tempprofile = CryoGrid.Presets.SamoylovDefault
grid = CryoGrid.Presets.DefaultGrid_5cm
initT = initializer(:T, tempprofile)
tile = CryoGrid.Presets.SoilHeatTile(
    :H,
    TemperatureGradient(forcings.Tair, NFactor(nf=Param(0.5), nt=Param(0.9))),
    GeothermalHeatFlux(0.053u"W/m^2"),
    soilprofile,
    initT;
    grid=grid
)
tspan = (DateTime(2010,10,30),DateTime(2011,10,30))
u0, du0 = initialcondition!(tile, tspan)
# Collect model parameters
p = CryoGrid.parameters(tile)
# Set up the `CryoGridProblem`
prob = CryoGridProblem(tile, u0, tspan, p, savevars=(:T,), step_limiter=nothing)

function loss(p)
    u0, du0 = initialcondition!(tile, tspan, p)
    prob = CryoGridProblem(tile, u0, tspan, p, savevars=(:T,), step_limiter=nothing)
    sol = @time solve(prob, Euler(), dt=300.0, saveat=24*3600.0, progress=true);
    return mean(sol.u[end])
end

using ForwardDiff
pvec = collect(p)
grad = ForwardDiff.gradient(loss, pvec)
@show grad
