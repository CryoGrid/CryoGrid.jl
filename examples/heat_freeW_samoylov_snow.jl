using CryoGrid
using CryoGrid.Physics
using Plots

forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044, :Tair => u"°C", :snowfall => u"mm/d");
# use air temperature as upper boundary forcing;
tair = TimeSeriesForcing(forcings.data.Tair, forcings.timestamps, :Tair);
# swe = TimeSeriesForcing(ustrip.(forcings.data.swe), forcings.timestamps, :swe);
# ρsn = TimeSeriesForcing(ustrip.(forcings.data.ρsn), forcings.timestamps, :ρsn);
snowfall = TimeSeriesForcing(uconvert.(u"m/s", forcings.data.snowfall.*1.0), forcings.timestamps, :snowfall)
# use default profiles for samoylov
soilprofile, tempprofile = CryoGrid.Presets.SamoylovDefault
# "simple" heat conduction model w/ 5 cm grid spacing (defaults to free water freezing scheme)
modelgrid = Grid(vcat([-1.0u"m"], CryoGrid.Presets.DefaultGrid_5cm))
initT = initializer(:T, tempprofile)
z_top = -2.0u"m"
z_sub = map(knot -> knot.depth, soilprofile)
z_bot = modelgrid[end]
snowmass = SnowMassBalance(
    para = Snow.DynamicSnow(
        ablation = Snow.DegreeDayMelt(factor=5.0u"mm/K/d")
    )
)
strat = @Stratigraphy(
    z_top => Top(TemperatureGradient(tair), Snowfall(snowfall)),
    z_top => :snowpack => Snowpack(Coupled(snowmass, HeatBalance()), para=Snow.Bulk(thresh=2.0u"cm")),
    z_sub[1] => :topsoil1 => Soil(HeatBalance(), para=soilprofile[1].value),
    z_sub[2] => :topsoil2 => Soil(HeatBalance(), para=soilprofile[2].value),
    z_sub[3] => :sediment1 => Soil(HeatBalance(), para=soilprofile[3].value),
    z_sub[4] => :sediment2 => Soil(HeatBalance(), para=soilprofile[4].value),
    z_sub[5] => :sediment3 => Soil(HeatBalance(), para=soilprofile[5].value),
    z_bot => Bottom(GeothermalHeatFlux(0.053u"J/s/m^2"))
);
tile = Tile(strat, modelgrid, initT)
# define time span, 2 years + 3 months
tspan = (DateTime(2016,9,30),DateTime(2018,12,31))
u0, du0 = initialcondition!(tile, tspan)
prob = CryoGridProblem(tile, u0, tspan, saveat=24*3600.0, savevars=(:T,:snowpack => (:dsn,:T_ub)))
# forward Euler with initial timestep of 5 minutes
sol = @time solve(prob, Euler(), dt=300.0, saveat=24*3600.0, progress=true);
out = CryoGridOutput(sol)
# Plot it!
zs = [1,10,20,30,50,100,200,500,1000]u"cm"
cg = Plots.cgrad(:copper,rev=true);
plot(ustrip(out.T[Z(Near(zs))]), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature (°C)", leg=false, dpi=150)
plt1 = plot!(ustrip(out.snowpack.T_ub), color=:skyblue, linestyle=:dash, alpha=0.7, leg=false, dpi=150)
plot(ustrip(out.swe), ylabel="Depth (m)", label="Snow water equivalent", dpi=150)
plt2 = plot!(ustrip(out.snowpack.dsn), label="Snow depth", legendtitle="", dpi=150)
plot(plt1, plt2, size=(1200,400), margins=5*Plots.Measures.mm)