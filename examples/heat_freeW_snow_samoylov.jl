# # [Soil heat with bulk snow scheme](@id example4)
# In this example, we construct a `Tile` consisting of a soil column with (i) heat conduction
# and (ii) a bulk (single-layer) snow scheme. The snowfall data comes
# from the ERA-Interim reanalysis product.

using CryoGrid
using OrdinaryDiffEq

# First we set up the model:
forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044);
soilprofile, tempprofile = CryoGrid.Presets.SamoylovDefault
initT = initializer(:T, tempprofile)
initsat = initializer(:sat, 1.0)
z_top = -2.0u"m"
z_sub = keys(soilprofile)
z_bot = 1000.0u"m"
upperbc = WaterHeatBC(
    SurfaceWaterBalance(forcings),
    TemperatureBC(forcings.Tair)
)
snowmass = SnowMassBalance(
    ablation = Snow.DegreeDayMelt(factor=5.0u"mm/K/d")
)
snowpack = Snowpack(
    para=Snow.Bulk(thresh=2.0u"cm"),
    mass=snowmass,
)
ground_layers = map(soilprofile) do para
    Ground(para, heat=HeatBalance(), water=WaterBalance())
end
strat = @Stratigraphy(
    z_top => Top(upperbc),
    z_top => snowpack,
    ground_layers...,
    z_bot => Bottom(GeothermalHeatFlux(0.053u"J/s/m^2"))
);
modelgrid = CryoGrid.Presets.DefaultGrid_5cm
tile = Tile(strat, modelgrid, initT, initsat)
# define time span, 2 years + 3 months
tspan = (DateTime(2010,9,30), DateTime(2012,9,30))
u0, du0 = @time initialcondition!(tile, tspan)
prob = CryoGridProblem(tile, u0, tspan, saveat=3*3600.0, savevars=(:T,:snowpack => (:dsn,:T_ub)))

# set up integrator
integrator = init(prob, Euler(), dt=300.0)
# advance 24 hours for testing
@time step!(integrator, 24*3600.0)

state = getstate(integrator)
state.snowpack.θwi

# solve full tspan with forward Euler and initial timestep of 5 minutes
@info "Running model ..."
sol = @time solve(prob, Euler(), dt=300.0, saveat=3*3600.0);
out = CryoGridOutput(sol)

# Plot it!
using Plots: plot, plot!, heatmap, cgrad, Measures
zs = [1,10,20,30,50,100,200,500]u"cm"
cg = cgrad(:copper,rev=true);
plot(ustrip(out.T[Z(Near(zs))]), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature (°C)", leg=false, dpi=150)
plt1 = plot!(ustrip.(out.snowpack.T_ub), color=:skyblue, linestyle=:dash, alpha=0.7, leg=false, dpi=150)

# Plot snow water equivalent and depth:
plot(ustrip(out.snowpack.swe), ylabel="Depth (m)", label="Snow water equivalent", dpi=150)
plt2 = plot!(ustrip.(out.snowpack.dsn), label="Snow depth", legend=nothing, legendtitle=nothing, dpi=150)
plot(plt1, plt2, size=(1600,700), margins=5*Measures.mm)

# Temperature heatmap:
T_sub = out.T[Z(Between(0.0u"m",10.0u"m"))]
heatmap(T_sub, yflip=true, size=(1200,600), dpi=150)

# Thaw depth:
td = Diagnostics.thawdepth(out.T)
plot(td, yflip=true, ylabel="Thaw depth (m)", size=(1200,600))

# ...and finally active layer thickness
alt = Diagnostics.active_layer_thickness(out.T)
plot(ustrip.(alt), ylabel="Active layer thickness (m)", xlabel="Number of years", label="ALT", size=(1200,600))
