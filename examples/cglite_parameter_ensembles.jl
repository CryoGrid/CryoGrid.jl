# # [Parameter ensembles](@id example12)
# This example is very similar to [Example 7](@ref) but shows how to create and run parallelized parameter ensembels.

# Make sure to explicitly import the `LiteImplicit` submodule which has
# the relevant solver types.
using CryoGrid
using CryoGrid.LiteImplicit
using Distributions

import Plots
import Random

if Threads.nthreads() == 1
    @warn "Only one thread is available. Ensemble execution will run sequentially. Did you start julia with `--threads=auto` ?"
end

# Load forcings and build stratigraphy like before, except this time we assign
# `Param` values to the quantiies which we want to vary in the ensemble. Here
# we vary the porosity in each layer as well as the n-factors.
forcings = loadforcings(CryoGrid.Forcings.Samoylov_ERA_MkL3_CCSM4_long_term);
soilprofile = SoilProfile(
    0.0u"m" => SimpleSoil(por=Param(0.80, prior=Uniform(0.65,0.95)),sat=1.0,org=0.75),
    0.1u"m" => SimpleSoil(por=Param(0.80, prior=Uniform(0.65,0.95)),sat=1.0,org=0.25),
    0.4u"m" => SimpleSoil(por=Param(0.55, prior=Uniform(0.35,0.75)),sat=1.0,org=0.25),
    3.0u"m" => SimpleSoil(por=Param(0.50, prior=Uniform(0.30,0.70)),sat=1.0,org=0.0),
    10.0u"m" => SimpleSoil(por=Param(0.30, prior=Uniform(0.10,0.50)),sat=1.0,org=0.0),
    100.0u"m" => SimpleSoil(por=0.05,sat=1.0,org=0.0),
)
z_top = -2.0u"m"
z_bot = 1000.0u"m"
upperbc = TemperatureBC(
    Input(:Tair),
    NFactor(
        nf=Param(0.5, prior=Beta(1,1)),
        nt=Param(0.9, prior=Beta(1,1)),
    ),
)
ssinit = ThermalSteadyStateInit(T0=-15.0u"°C")
heatop = Heat.EnthalpyImplicit()
heat = HeatBalance(heatop)
soil_layers = map(para -> Ground(para; heat), soilprofile)
strat = Stratigraphy(
    z_top => Top(upperbc),
    soil_layers,
    z_bot => Bottom(GeothermalHeatFlux(0.053u"W/m^2"))
);
modelgrid = CryoGrid.DefaultGrid_2cm
tile = Tile(strat, modelgrid, forcings, ssinit);
# Since the solver can take daily timesteps, we can easily specify longer simulation time spans at minimal cost.
# Here we specify a time span of 10 years.
tspan = (DateTime(2000,1,1), DateTime(2010,12,31))
u0, du0 = initialcondition!(tile, tspan);
prob = CryoGridProblem(tile, u0, tspan, saveat=24*3600.0, savevars=(:T,))

# Here we retrieve the `CryoGridParams` from the `CryoGridProblem` constructed above.
# The CryoGridParams type behaves like a table and can be easily converted
# to a DataFrame with DataFrame(params) when DataFrames.jl is loaded.
params = prob.p

# Note that you can also use Julia's `vec` method to convert `CryoGridParams` into a `ComponentVector` with labels.
p0 = vec(params)

# Here we extract prior distributions and collect them into a multivariate Product distribution;
# note that this assumes each parameter to be independent from the others
prior = Product(collect(params[:prior]))

# Method 1: SciML EnsembleProblem

function make_prob_func(ensmeble::AbstractMatrix)
    function prob_func(prob, i, repeat)
        return remake(prob, p=ensmeble[:,i])
    end
end

function output_func(sol, i)
    # return CryoGridOutput; false indicates that the run does not need to be repeated
    return CryoGridOutput(sol), false
end

# Now we sample parameter values from prior with fixed RNG;
# the number of samples determines the size of the ensemble
const rng = Random.MersenneTwister(1234)
prior_ensemble = rand(rng, prior, 64)

# We create an `EnsembleProblem` from `CryoGridProblem` and prob/output functions;
# note that we use safetycopy=true because we're using multithreading;
# this prevents the different threads from using the same state caches
prob_func = make_prob_func(prior_ensemble)
ensprob = EnsembleProblem(prob; prob_func, output_func, safetycopy=true)
# solve each trajectory with LiteImplicitEuler and multithreading;
# alternatively, one can specify EnsembleDistributed() for process or slurm parallelization or EnsembleSerial() for sequential execution.
enssol = @time solve(ensprob, LiteImplicitEuler(), EnsembleThreads(), trajectories=size(prior_ensemble,2))

# Now we will extract permafrost temperatures at 20m depth and plot the ensemble.
T20m_ens = reduce(hcat, map(out -> out.T[Z(Near(20.0u"m"))], enssol))
Plots.plot(T20m_ens, leg=nothing, c=:black, alpha=0.5, ylabel="Permafrost temperature")

alt_ens = reduce(hcat, map(out -> Diagnostics.active_layer_thickness(out.T), enssol))
Plots.plot(alt_ens, leg=nothing, c=:black, alpha=0.5, ylabel="Active layer thickness")

# Method 2: Simple (sequential) for-loop; this could also be parallelized with pmap or @distributed
# multithreading with @threads is also possible, but in this case, one should call deepcopy(prob) to prevent cache collisions

using ProgressMeter

@showprogress for i in axes(prior_ensemble)[2]
    p_i = prior_ensemble[:,i]
    # solve with parameters p_i
    sol_i = solve(prob, LiteImplicitEuler(), p=p_i)
    # here we would need to process the output and store the results;
    # omitted in this example for brevity
end
