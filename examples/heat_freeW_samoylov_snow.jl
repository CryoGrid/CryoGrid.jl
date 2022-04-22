using CryoGrid
using Plots

forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA5_fitted_daily_1979_2020, :Tair => u"°C", :Dsn => u"m"; spec=JsonSpec{2});
# use air temperature as upper boundary forcing;
tair = TimeSeriesForcing(ustrip.(forcings.data.Tair), forcings.timestamps, :Tair);
snowdepth = TimeSeriesForcing(ustrip.(forcings.data.Dsn), forcings.timestamps, :Dsn);
# use default profiles for samoylov
soilprofile, tempprofile = CryoGrid.Presets.SamoylovDefault
# "simple" heat conduction model w/ 5 cm grid spacing (defaults to free water freezing scheme)
modelgrid = Grid(vcat([-1.0u"m"], CryoGrid.Presets.DefaultGrid_5cm))
initT = initializer(:T, tempprofile)
z_top = -2.0u"m"
z_sub = map(knot -> knot.depth, soilprofile)
z_bot = modelgrid[end]
strat = @Stratigraphy(
    z_top => top(TemperatureGradient(tair)),
    z_top => subsurface(:snowpack, Snowpack(para=Snow.Bulk(dsn=snowdepth)), Heat(:H)),
    z_sub[1] => subsurface(:topsoil1, Soil(para=soilprofile[1].value), Heat(:H)),
    z_sub[2] => subsurface(:topsoil2, Soil(para=soilprofile[2].value), Heat(:H)),
    z_sub[3] => subsurface(:sediment1, Soil(para=soilprofile[3].value), Heat(:H)),
    z_sub[4] => subsurface(:sediment2, Soil(para=soilprofile[4].value), Heat(:H)),
    z_sub[5] => subsurface(:sediment3, Soil(para=soilprofile[5].value), Heat(:H)),
    z_bot => bottom(GeothermalHeatFlux(0.053u"J/s/m^2"))
);
tile = Tile(strat, modelgrid, initT)
# define time span
tspan = (DateTime(2010,10,30),DateTime(2012,12,30))
p = parameters(tile)
u0, du0 = initialcondition!(tile, tspan, p)
# CryoGrid front-end for ODEProblem
prob = CryoGridProblem(tile,u0,tspan,p,savevars=(:T,:C,:kc))
# solve with forward Euler, 10-minute time steps
out = @time solve(prob, Euler(), dt=600.0, saveat=24*3600.0, progress=true) |> CryoGridOutput;
# Plot it!
zs = [1,10,20,30,50,100,200,500,1000]u"cm"
cg = Plots.cgrad(:copper,rev=true);
plot(ustrip(out.T[Z(Near(zs))]), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", leg=false, dpi=150)
plot(ustrip(out.H[Z(Near(zs))]), color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Enthalpy", leg=false, dpi=150)

integrator = init(prob, Euler(), dt=300.0)
step!(integrator)
snow_state = getstate(:snowpack, integrator)
snow_state.kc
