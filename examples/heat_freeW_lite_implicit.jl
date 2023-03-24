using CryoGrid
using CryoGrid.LiteImplicit
using Plots

forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_MkL3_CCSM4_long_term, :Tair => u"°C", :Ptot => u"mm");
# use air temperature as upper boundary forcing 
tair = TimeSeriesForcing(forcings.data.Tair, forcings.timestamps, :Tair);
tempprofile_linear = TemperatureProfile(
    0.0u"m" => -30.0u"°C",
    10.0u"m" => -10.0u"°C", 
    1000.0u"m" => 10.2u"°C"
)
z_top = -2.0u"m"
z_bot = 1000.0u"m"
upperbc = TemperatureGradient(tair, NFactor())
initT = initializer(:T, tempprofile_linear)
@info "Building stratigraphy"
heatop = Heat.EnthalpyImplicit()
strat = @Stratigraphy(
    z_top => Top(upperbc),
    0.0u"m" => :topsoil1 => Soil(HomogeneousMixture(por=0.80,sat=1.0,org=0.75), heat=HeatBalance(heatop)),
    0.1u"m" => :topsoil2 => Soil(HomogeneousMixture(por=0.80,sat=1.0,org=0.25), heat=HeatBalance(heatop)),
    0.4u"m" => :sediment1 => Soil(HomogeneousMixture(por=0.55,sat=1.0,org=0.25), heat=HeatBalance(heatop)),
    3.0u"m" => :sediment2 => Soil(HomogeneousMixture(por=0.50,sat=1.0,org=0.0), heat=HeatBalance(heatop)),
    10.0u"m" => :sediment3 => Soil(HomogeneousMixture(por=0.30,sat=1.0,org=0.0), heat=HeatBalance(heatop)),
    z_bot => Bottom(GeothermalHeatFlux(0.053u"W/m^2"))
);
@info "Building tile"
tile = @time Tile(strat, AutoGrid(), initT)
# define time span, 5 years
tspan = (DateTime(2010,12,30), DateTime(2015,12,30))
tspan_sol = convert_tspan(tspan)
u0, du0 = @time initialcondition!(tile, tspan);
prob = CryoGridProblem(tile, u0, tspan, saveat=24*3600.0, savevars=(:θw,:T,))
@info "Running model"
sol = @time solve(prob, LiteImplicitEuler(), dt=24*3600)
out = CryoGridOutput(sol)
# Plot the results
zs = [5,10,15,20,25,30,40,50,100,500,1000,5000]u"cm"
cg = Plots.cgrad(:copper,rev=true);
plot(out.T[Z(Near(zs))] |> ustrip, color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", title="", leg=false, dpi=150)
