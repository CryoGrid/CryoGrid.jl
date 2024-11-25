using CryoGrid
using CryoGrid.Presets
using Dates
using Plots
using Statistics
using BenchmarkTools
using Serialization

function run_solver_benchmark(solvers, grids, dts, tspan, upperbc;
    freezecurve=FreeWater(), params=[], dt_ref=30.0, abstol=1e-3, reltol=1e-4, adaptive=true, savefreq=6*3600.0, profile=SamoylovDefault)
    errors = zeros(length(grids), length(solvers), 2)
    times = zeros(length(grids), length(solvers))
    for (i, (dt, grid)) in enumerate(zip(dts,grids))
        min_dx = round(minimum(Δ(grid)), digits=2)
        # high accuracy reference solution, 1 minute time steps
        println("Computing reference solution for grid spacing $min_dx with forward Euler @ $dt_ref s time steps ...")
        # basic 1-layer heat conduction model (defaults to free water freezing scheme)
        model = Presets.SoilHeat(:H, upperbc, profile; grid=grid, freezecurve=freezecurve)
        p = copy(model.pproto)
        p .= params
        prob = CryoGridProblem(model,tspan,p)
        out_ref = solve(prob, Euler(), dt=30.0, saveat=6*3600.0, progress=true) |> CryoGridOutput
        # max_alt = findlast(any(out_ref.soil.T >= 0.0, dims=2))
        # err_lb = max_alt + 10 # set lower bound for error to be max ALT + 10 grid cells
        for (j, solver) in enumerate(solvers)
            println("$solver w/ minimum grid spacing $(min_dx)")
            # warm-up
            println("Running warm-up ...")
            prob = CryoGridProblem(model,(tspan[1], tspan[1] + Day(1)), p)
            integrator = init(prob, solver(), dt=dt, abstol=abstol, reltol=reltol, adaptive=adaptive, saveat=6*3600.0, progress=true)
            @time step!(integrator, prob.tspan[end] - prob.tspan[1])
            prob = remake(prob, tspan=Dates.datetime2epochms.(tspan)./1000.0)
            println("Running full time span ...")
            integrator = init(prob, solver(), dt=dt, abstol=abstol, reltol=reltol, adaptive=adaptive, saveat=6*3600.0, progress=true)
            Δt = prob.tspan[end] - prob.tspan[1]
            res = @benchmark step!($integrator, $Δt) samples=1
            @show res
            out = CryoGridOutput(integrator.sol)
            # absolute error integrated over space
            R = abs.(out.soil.H .- out_ref.soil.H)
            err = R'*Δ(grid)
            # compute mean + stderr of integrated absolute error over time, ignore first 24 hours
            errors[i, j, 1] = err_mean = mean(err[4:end])
            errors[i, j, 2] = err_std = std(err[4:end])
            times[i, j] = res.times[1] / 1e9 # in seconds
            println("Done! H error: $err_mean +/- $err_std")
        end
    end
    return (runtime=times, error=errors)
end

# Test 1: Energy + air temperature upper bc + free water freeze curve
filename = "./input/FORCING_JSONfiles/samoylov_ERA_obs_fitted_1979_2014_spinup_extended2044.json"
forcings = loadforcings(filename);
# use air temperature as upper boundary forcing
tair = TemperatureBC(forcings.Tair)
solvers = [Euler, DP5, ROCK2, ROCK4, Trapezoid, ROS3P]
grids = [Presets.DefaultGrid_2cm, Presets.DefaultGrid_5cm, Presets.DefaultGrid_10cm, Presets.DefaultGrid_20cm]
dts = [2*60.0, 10*60.0, 30*60.0, 3600.0]
tspan = (DateTime(2010,10,30),DateTime(2011,10,30))
results_freeW_tair = run_solver_benchmark(solvers, grids, dts, tspan, tair)
serialize("solver_benchmark_freeW_tair.ser", results_freeW_tair)

results_freeW_tair = deserialize("solver_benchmark_freeW_tair.ser")
scatter([2, 5, 10, 20], results_freeW_tair.runtime, labels=["Euler" "DP5" "ROCK2" "ROCK4" "Crank-Nicolson" "ROS3P"], xlabel="Minimum grid spacing (cm)",
    ylabel="Run time per simulation year (s)", title="Integrator run times w/ free water FC, air temp. bc", dpi=150)
hline!([1.0], linestyle=:dash, c=:black, label=nothing)

xnames = repeat(["0.02", "0.05", "0.10", "0.20"], outer=length(solvers))
groupnames = repeat([string(s) for s in solvers], inner=length(grids))
StatsPlots.groupedbar(xnames, results_freeW_tair.runtime, group=groupnames, 
    bar_position=:dodge, ylabel="Run time per simulation year (s)", xlabel="Minimum grid spacing (m)",
    title="Run times, free water FC, air temp. upper bc")
savefig("solver_benchmark_freeW_tair_runtimes.png")
StatsPlots.groupedbar(xnames, results_freeW_tair.error[:,:,1], yerror=results_freeW_tair.error[:,:,2], group=groupnames, 
    bar_position=:dodge, ylims=(0,6e5), ylabel="Error w.r.t reference run (J/m³)", xlabel="Minimum grid spacing (m)",
    title="Error, free water FC, air temp. upper bc", legend=:topleft)
savefig("solver_benchmark_freeW_tair_error.png")

# Test 2: Energy + SEB upper bc + free water freeze curve
forcings = (; Tair,pr,q,wind,Lin,Sin)
seb = SurfaceEnergyBalance()
solvers = [Euler, DP5, ROCK2, ROCK4]
results_freeW_seb = run_solver_benchmark(solvers, grids, dts, tspan, seb)
serialize("solver_benchmark_freeW_seb.ser", results_freeW_seb)

results_freeW_seb = deserialize("solver_benchmark_freeW_seb.ser")
scatter([2, 5, 10, 20], results_freeW_seb.runtime, labels=["Euler" "DP5" "ROCK2" "ROCK4"], xlabel="Minimum grid spacing (cm)",
    ylabel="Run time per simulation year (s)", title="Integrator run times w/ free water FC, SEB bc", dpi=150)
hline!([1.0], linestyle=:dash, c=:black, label=nothing)

plot([2, 5, 10, 20], results_freeW_seb.error[:,:,1], ribbon=results_freeW_seb.error[:,:,2].*sqrt(4*365),
    labels=["Euler" "DP5" "ROCK2" "ROCK4"], xlabel="Minimum grid spacing (cm)", ylabel="Error (J/m³)",
    title="Integrator error w/ free water FC + SEB bc", legend=:topleft, w=3, fillalpha=0.2, dpi=100)

xnames = repeat(["0.02", "0.05", "0.10", "0.20"], outer=length(solvers))
groupnames = repeat([string(s) for s in solvers], inner=length(grids))
StatsPlots.groupedbar(xnames, results_freeW_seb.runtime, group=groupnames, 
    bar_position=:dodge, ylabel="Run time per simulation year (s)", xlabel="Minimum grid spacing (m)",
    title="Run times, free water FC, SEB upper bc")
savefig("solver_benchmark_freeW_seb_runtimes.png")
StatsPlots.groupedbar(xnames, results_freeW_seb.error[:,:,1], yerror=results_freeW_seb.error[:,:,2], group=groupnames, 
    bar_position=:dodge, ylims=(0,6e5), ylabel="Error w.r.t reference run (J/m³)", xlabel="Minimum grid spacing (m)",
    title="Error, free water FC, SEB upper bc", legend=:topleft)
savefig("solver_benchmark_freeW_seb_error.png")

# Test 3: Energy + SEB upper bc + Dall'Amico freeze curve
solvers = [Euler]
p = [4.0, 2.0, 0.0]
results_vgfc_seb = run_solver_benchmark(solvers, grids, dts, tspan, seb, params=p, adaptive=false, freezecurve=SFCC(DallAmico()))
serialize("solver_benchmark_vgfc_seb.ser", results_vgfc_seb)

results_vgfc_seb = deserialize("solver_benchmark_vgfc_seb.ser")
plot([2, 5, 10, 20], results_vgfc_seb.runtime, labels="Euler", xlabel="Minimum grid spacing (cm)",
    ylabel="Run time per simulation year (s)", title="Integrator run times w/ Dall'Amico FC + SEB bc", dpi=150)
hline!([1.0], linestyle=:dash, c=:black, label=nothing)

xnames = repeat(["0.02", "0.05", "0.10", "0.20"], outer=length(solvers))
groupnames = repeat([string(s) for s in solvers], inner=length(grids))
StatsPlots.groupedbar(xnames, results_vgfc_seb.runtime, group=groupnames, 
    bar_position=:dodge, ylabel="Run time per simulation year (s)", xlabel="Minimum grid spacing (m)",
    title="Run times, Dall'Amico FC, SEB upper bc")
savefig("solver_benchmark_vgc_seb_runtimes.png")
StatsPlots.groupedbar(xnames, results_vgfc_seb.error[:,:,1], yerror=results_vgfc_seb.error[:,:,2], group=groupnames, 
    bar_position=:dodge, ylims=(0,6e5), ylabel="Error w.r.t reference run (J/m³)", xlabel="Minimum grid spacing (m)",
    title="Error, Dall'Amico FC, SEB upper bc", legend=:topleft)
savefig("solver_benchmark_vgfc_seb_error.png")
