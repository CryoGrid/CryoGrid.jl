using CryoGrid
using Plots
using Statistics

forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044, :Tair => u"°C");
# use air temperature as upper boundary forcing;
tair = TimeSeriesForcing(forcings.data.Tair, forcings.timestamps, :Tair);
# use default profiles for samoylov
soilprofile, tempprofile = CryoGrid.Presets.SamoylovDefault
# "simple" heat conduction model w/ 5 cm grid spacing (defaults to free water freezing scheme)
grid = CryoGrid.Presets.DefaultGrid_5cm
initT = initializer(:T, tempprofile)
# construct soil Tile with heat conduction and parameterized n-factors
tile = CryoGrid.Presets.SoilHeatTile(
    :H,
    TemperatureGradient(tair, NFactor(nf=Param(0.5), nt=Param(0.9))),
    GeothermalHeatFlux(0.053u"W/m^2"),
    soilprofile,
    initT;
    grid=grid
)
# define time span
tspan = (DateTime(2010,10,30),DateTime(2011,10,30))
u0, du0 = initialcondition!(tile, tspan)
p = CryoGrid.parameters(tile)
# construct CryoGridProblem with tile, initial condition, and timespan;
# we disable the default timestep limiter since we will use an adaptive solver.
prob = CryoGridProblem(tile, u0, tspan, collect(p), savevars=(:T,), step_limiter=nothing)

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
