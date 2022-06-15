using CryoGrid
using Plots

# forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA5_fitted_daily_1979_2020, :Tair => u"°C", :swe => u"m", :ρsn => u"kg/m^3"; spec=JsonSpec{2});
forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044, :Tair => u"°C", :snowfall => u"mm/d"; spec=JsonSpec{1});
# use air temperature as upper boundary forcing;
tair = TimeSeriesForcing(ustrip.(forcings.data.Tair), forcings.timestamps, :Tair);
# swe = TimeSeriesForcing(ustrip.(forcings.data.swe), forcings.timestamps, :swe);
# ρsn = TimeSeriesForcing(ustrip.(forcings.data.ρsn), forcings.timestamps, :ρsn);
snowfall = TimeSeriesForcing(convert.(Float64, ustrip.(u"m/s", forcings.data.snowfall)), forcings.timestamps, :snowfall)
# use default profiles for samoylov
soilprofile, tempprofile = CryoGrid.Presets.SamoylovDefault
# "simple" heat conduction model w/ 5 cm grid spacing (defaults to free water freezing scheme)
modelgrid = Grid(vcat([-1.0u"m"], CryoGrid.Presets.DefaultGrid_5cm))
initT = initializer(:T, tempprofile)
z_top = -2.0u"m"
z_sub = map(knot -> knot.depth, soilprofile)
z_bot = modelgrid[end]
snowmass = SnowMassBalance(
    para = Snow.Dynamic(
        ablation = Snow.DegreeDayMelt(factor=5.0u"mm/K/d")
    )
)
strat = @Stratigraphy(
    z_top => top(TemperatureGradient(tair), Snowfall(snowfall)),
    # prescribed snow
    # z_top => subsurface(:snowpack, Snowpack(para=Snow.Bulk()), SnowMassBalance(para=Snow.Prescribed(swe=swe, ρsn=ρsn)), Heat(:H)),
    # "dynamic" snow (i.e. modeled snow accumulation and ablation)
    z_top => subsurface(:snowpack, Snowpack(para=Snow.Bulk(thresh=2.0u"cm")), snowmass, Heat(:H)),
    z_sub[1] => subsurface(:topsoil1, Soil(para=soilprofile[1].value), Heat(:H)),
    z_sub[2] => subsurface(:topsoil2, Soil(para=soilprofile[2].value), Heat(:H)),
    z_sub[3] => subsurface(:sediment1, Soil(para=soilprofile[3].value), Heat(:H)),
    z_sub[4] => subsurface(:sediment2, Soil(para=soilprofile[4].value), Heat(:H)),
    z_sub[5] => subsurface(:sediment3, Soil(para=soilprofile[5].value), Heat(:H)),
    z_bot => bottom(GeothermalHeatFlux(0.053u"J/s/m^2"))
);
tile = Tile(strat, modelgrid, initT)
# define time span
tspan = (DateTime(2010,9,30),DateTime(2012,9,30))
p = parameters(tile)
u0, du0 = initialcondition!(tile, tspan, p)
prob = CryoGridProblem(tile,u0,tspan,p,step_limiter=nothing,savevars=(:T,:snowpack => (:dsn,:T_ub)))
sol = @time solve(prob, SSPRK22(), dt=300.0, saveat=24*3600.0, progress=true);
out = CryoGridOutput(sol)
# Plot it!
zs = [1,10,20,30,50,100,200,500,1000]u"cm"
cg = Plots.cgrad(:copper,rev=true);
plot(ustrip(out.T[Z(Near(zs))]), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature (°C)", leg=false, dpi=150)
plt1 = plot!(ustrip(out.snowpack.T_ub), color=:skyblue, linestyle=:dash, alpha=0.7, leg=false, dpi=150)
plot(ustrip(out.swe), ylabel="Depth (m)", label="Snow water equivalent", dpi=150)
plt2 = plot!(ustrip(out.snowpack.dsn), label="Snow depth", legendtitle="", dpi=150)
plot(plt1, plt2, size=(1200,400), margins=5*Plots.Measures.mm)
