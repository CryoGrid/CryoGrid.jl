# # [Soil heat w/ SEB, snow cover, and bucket water scheme](@id example6)
# In this example, we construct a `Tile` consisting of a soil column with (i) heat conduction
# forced by the surface energy balance (SEB), (ii) a bulk snow scheme, and
# (iii) a bucket hydrology scheme.

# For this example, we need to use an OrdinaryDiffEq integrator.
using CryoGrid
using OrdinaryDiffEq

# First, load the forcings and construct the Tile.
modelgrid = CryoGrid.Presets.DefaultGrid_2cm;
soilprofile = SoilProfile(
    0.0u"m" => SimpleSoil(por=0.80,sat=1.0,org=0.75),
    0.1u"m" => SimpleSoil(por=0.80,sat=1.0,org=0.25),
    0.4u"m" => SimpleSoil(por=0.55,sat=1.0,org=0.25),
    3.0u"m" => SimpleSoil(por=0.50,sat=1.0,org=0.0),
    10.0u"m" => SimpleSoil(por=0.30,sat=1.0,org=0.0),
);
## mid-winter temperature profile
tempprofile = CryoGrid.Presets.SamoylovDefault.tempprofile
forcings = loadforcings(CryoGrid.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044);
tempprofile = CryoGrid.Presets.SamoylovDefault.tempprofile
initT = initializer(:T, tempprofile)
seb = SurfaceEnergyBalance()
swb = SurfaceWaterBalance()
upperbc = WaterHeatBC(swb, seb)
heat = HeatBalance()
water = WaterBalance(BucketScheme(), DampedET())
## build stratigraphy
strat = @Stratigraphy(
    -z => Top(upperbc), 
    -z => Snowpack(heat=HeatBalance(), water=water),
    0.0u"m" => Ground(soilprofile; heat, water),
    1000.0u"m" => Bottom(GeothermalHeatFlux(0.053u"J/s/m^2")),
);
## create Tile
tile = Tile(strat, modelgrid, forcings, initT);

# Set up the problem and solve it!
tspan = (DateTime(2010,10,30), DateTime(2011,10,30))
## generate initial condition and set up CryoGridProblem
u0, du0 = @time initialcondition!(tile, tspan)

prob = CryoGridProblem(
    tile,
    u0,
    tspan,
    savevars=(:T,:jH,:top => (:Qh,:Qe,:Qg,),:snowpack => (:dsn,)),
    saveat=3*3600.0
)
integrator = init(prob, Euler(), dt=60.0)
## step forwards 24 hours and check for NaN/Inf values
@time step!(integrator); integrator.dt
@assert all(isfinite.(integrator.u))
## iterate over remaining timespan at fixed points using `TimeChoiceIterator`
@time for (u,t) in TimeChoiceIterator(integrator, convert_t.(tspan[1]:Day(1):tspan[end]))
    state = getstate(integrator)
    @assert isfinite(state.top.Qg[1])
    @info "Current t=$(Date(convert_t(t))), dt=$(integrator.dt), Tsurf=$(state.ground.T[1]), infil=$(state.top.jw_infil[1])"
    if integrator.sol.retcode != ReturnCode.Default
        break
    end
end
out = CryoGridOutput(integrator.sol)

# Plot it!
import Plots
zs = [1,3,5,7,10,15,20,25,30,40,50,100,150,200,500,1000]u"cm"
cg = Plots.cgrad(:copper,rev=true);
Plots.plot(ustrip.(out.T[Z(Near(zs))]), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, size=(800,500), dpi=150)

# Saturation:
Plots.plot(ustrip.(out.sat[Z(Near(zs))]), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Soil saturation", leg=false, size=(800,500), dpi=150)

# Runoff
Plots.plot(ustrip.(out.top.runoff), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Total runoff", leg=false, size=(800,500), dpi=150)

# Evapotranspiration
Plots.plot(ustrip.(out.top.ET), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Total ET", leg=false, size=(800,500), dpi=150)

# Snow depth:
Plots.plot(ustrip.(out.snowpack.dsn), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Snow depth", leg=false, size=(800,500), dpi=150)

# Integrated ground heat flux:
Plots.plot(ustrip.(cumsum(out.top.Qg, dims=2)), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Integrated ground heat flux", leg=false, size=(800,500), dpi=150)

# Integratoed ground latent heat flux:
Plots.plot(ustrip.(cumsum(out.top.Qe, dims=2)), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Integrated ground heat flux", leg=false, size=(800,500), dpi=150)

using BenchmarkTools
using Interpolations

function interpolant(xs::AbstractVector, ts::AbstractVector)
    f = interpolate((ts,), xs, Gridded(Linear()))
    return f
end

function eval(xs::AbstractVector, ts::AbstractVector, t)
    f = interpolate((ts,), xs, Gridded(Linear()))
    return f(t)
end

const xdata = randn(1000)
const tdata = 1:1000
@btime eval(xdata, tdata, 10.0)

const itp = interpolant(xdata, tdata)
@btime itp(10.0)

x = randn(1000,2)
interpolate((tdata,1:2), x, (Gridded(Linear()), NoInterp()))