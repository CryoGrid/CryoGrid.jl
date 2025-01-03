# # [Soil heat with bulk snow scheme](@id example4)
# In this example, we construct a `Tile` consisting of a soil column with (i) heat conduction
# and (ii) a bulk (single-layer) snow scheme. The snowfall data comes
# from the ERA-Interim reanalysis product.

using CryoGrid
using OrdinaryDiffEq

# First we set up the model:
forcings = loadforcings(CryoGrid.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044);
soilprofile = SoilProfile(
    0.0u"m" => SimpleSoil(por=0.80,sat=1.0,org=0.75),
    0.1u"m" => SimpleSoil(por=0.80,sat=1.0,org=0.25),
    0.4u"m" => SimpleSoil(por=0.55,sat=1.0,org=0.25),
    3.0u"m" => SimpleSoil(por=0.50,sat=1.0,org=0.0),
    10.0u"m" => SimpleSoil(por=0.30,sat=1.0,org=0.0),
);
initT = initializer(:T, CryoGrid.SamoylovDefault.tempprofile)
initsat = initializer(:sat, 1.0)
z_top = -2.0u"m"
z_sub = keys(soilprofile)
z_bot = 1000.0u"m"
upperbc = WaterHeatBC(
    SurfaceWaterBalance(),
    TemperatureBC(Input(:Tair))
)
snowpack = Snowpack(
    para=Snow.Bulk(),
    mass=SnowMassBalance(ablation = Snow.DegreeDayMelt(factor=5.0u"mm/K/d")),
    heat=HeatBalance(),
    water=WaterBalance(),
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
modelgrid = CryoGrid.DefaultGrid_5cm
tile = Tile(strat, modelgrid, forcings, initT, initsat)
# define time span, 2 years + 3 months
tspan = (DateTime(2010,9,30), DateTime(2012,9,30))
u0, du0 = @time initialcondition!(tile, tspan)
prob = CryoGridProblem(tile, u0, tspan, saveat=3*3600.0, savevars=(:T, :top => (:T_ub), :snowpack => (:dsn,)))

# solve full tspan with forward Euler and initial timestep of 5 minutes
@info "Running model ..."
sol = @time solve(prob, Euler(), dt=300.0);
out = CryoGridOutput(sol)

# Plot it!
using Plots: plot, plot!, heatmap, cgrad, Measures
zs = [1,10,20,30,50,100,200,500]u"cm"
cg = cgrad(:copper,rev=true);
plot(ustrip(out.T[Z(Near(zs))]), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature (°C)", leg=false, dpi=150)
plot!(ustrip(out.T[1,:]), color=:darkgray, ylabel="Temperature (°C)", leg=false, dpi=150)
plt1 = plot!(ustrip.(out.top.T_ub), color=:skyblue, linestyle=:dash, alpha=0.5, leg=false, dpi=150)

# Plot snow water equivalent and depth:
plot(ustrip(out.snowpack.swe), ylabel="Depth (m)", label="Snow water equivalent", dpi=150)
plt2 = plot!(ustrip.(out.snowpack.dsn), label="Snow depth", ylabel="Depth (m)", legendtitle=nothing, dpi=150)
plot(plt1, plt2, size=(1600,700), margins=5*Measures.mm)

# Temperature heatmap:
T_sub = out.T[Z(Between(0.0u"m",10.0u"m"))]
heatmap(T_sub, yflip=true, size=(1200,600), dpi=150)

# Thaw depth:
td = Diagnostics.thawdepth(out.T[Z(Where(>=(0.0u"m")))])
plot(td, yflip=true, ylabel="Thaw depth (m)", size=(1200,600))
