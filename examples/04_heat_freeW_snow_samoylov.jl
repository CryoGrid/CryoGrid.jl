# # Example 4
# ## Heat conduction on soil column with bulk snow scheme
# In this example, we construct a `Tile` consisting of a soil column with (i) heat conduction
# and (ii) a bulk (single-layer) snow scheme. The snowfall data comes
# from the ERA-Interim reanalysis product.

using CryoGrid

forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044);
soilprofile, tempprofile = CryoGrid.Presets.SamoylovDefault
initT = initializer(:T, tempprofile)
z_top = -2.0u"m"
z_sub = map(knot -> knot.depth, soilprofile)
z_bot = 1000.0u"m"
snowmass = SnowMassBalance(
    para = Snow.DynamicSnow(
        ablation = Snow.DegreeDayMelt(factor=5.0u"mm/K/d")
    )
)
upperbc = WaterHeatBC(
    SurfaceWaterBalance(rainfall=forcings.rainfall, snowfall=forcings.snowfall),
    TemperatureGradient(forcings.Tair)
)
strat = @Stratigraphy(
    z_top => Top(upperbc),
    z_top => :snowpack => Snowpack(para=Snow.Bulk(thresh=2.0u"cm"), mass=snowmass, heat=HeatBalance()),
    z_sub[1] => :topsoil1 => SimpleSoil(soilprofile[1].value, heat=HeatBalance()),
    z_sub[2] => :topsoil2 => SimpleSoil(soilprofile[2].value, heat=HeatBalance()),
    z_sub[3] => :sediment1 => SimpleSoil(soilprofile[3].value, heat=HeatBalance()),
    z_sub[4] => :sediment2 => SimpleSoil(soilprofile[4].value, heat=HeatBalance()),
    z_sub[5] => :sediment3 => SimpleSoil(soilprofile[5].value, heat=HeatBalance()),
    z_bot => Bottom(GeothermalHeatFlux(0.053u"J/s/m^2"))
);
modelgrid = CryoGrid.Presets.DefaultGrid_5cm
tile = Tile(strat, modelgrid, initT)
# define time span, 2 years + 3 months
tspan = (DateTime(2010,9,30), DateTime(2012,9,30))
u0, du0 = initialcondition!(tile, tspan)
prob = CryoGridProblem(tile, u0, tspan, saveat=3*3600.0, savevars=(:T,:snowpack => (:dsn,:T_ub)))

# set up integrator
integrator = init(prob, Euler(), dt=300.0, saveat=3*3600.0)
# advance 24 hours
@time step!(integrator, 24*3600.0)

# solve full tspan with forward Euler and initial timestep of 5 minutes
sol = @time solve(prob, Euler(), dt=300.0, saveat=3*3600.0, progress=true);
out = CryoGridOutput(sol)

# Plot it!
import Plots

zs = [1,10,20,30,50,100,200,500]u"cm"
cg = Plots.cgrad(:copper,rev=true);
Plots.plot(ustrip(out.T[Z(Near(zs))]), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature (°C)", leg=false, dpi=150)
plt1 = Plots.plot!(ustrip.(out.snowpack.T_ub), color=:skyblue, linestyle=:dash, alpha=0.7, leg=false, dpi=150)
Plots.plot(ustrip(out.snowpack.swe), ylabel="Depth (m)", label="Snow water equivalent", dpi=150)
plt2 = Plots.plot!(ustrip.(out.snowpack.dsn), label="Snow depth", legend=nothing, legendtitle=nothing, dpi=150)
Plots.plot(plt1, plt2, size=(1600,700), margins=5*Plots.Measures.mm)
# heatmap
T_sub = out.T[Z(Between(0.0u"m",10.0u"m"))]
Plots.heatmap(T_sub, yflip=true, size=(1200,600), dpi=150)
# thaw depth
td = Diagnostics.thawdepth(out.T)
Plots.plot(td, yflip=true, ylabel="Thaw depth (m)", size=(1200,600))
# active layer thickness
alt = Diagnostics.active_layer_thickness(out.T)
Plots.plot(ustrip.(alt), ylabel="Active layer thickness (m)", xlabel="Number of years", label="ALT", size=(1200,600))
