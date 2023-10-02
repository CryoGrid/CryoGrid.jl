using CryoGrid
using CryoGrid.LiteImplicit
using Plots

CryoGrid.debug(true)

forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_MkL3_CCSM4_long_term, :Tair => u"°C", :Ptot => u"mm");
# use air temperature as upper boundary forcing 
tair = TimeSeriesForcing(forcings.data.Tair, forcings.timestamps, :Tair);
# snowdepth = TimeSeriesForcing(ustrip.(forcings.data.Dsn), forcings.timestamps, :Dsn);
soilprofile = SoilProfile(
    0.0u"m" => HomogeneousMixture(por=0.80,sat=0.9,org=0.75), #(θwi=0.80,θm=0.05,θo=0.15,ϕ=0.80),
    0.1u"m" => HomogeneousMixture(por=0.80,sat=0.9,org=0.25), #(θwi=0.80,θm=0.15,θo=0.05,ϕ=0.80),
    0.4u"m" => HomogeneousMixture(por=0.55,sat=0.9,org=0.25), #(θwi=0.80,θm=0.15,θo=0.05,ϕ=0.55),
    3.0u"m" => HomogeneousMixture(por=0.50,sat=1.0,org=0.0), #(θwi=0.50,θm=0.50,θo=0.0,ϕ=0.50),
    10.0u"m" => HomogeneousMixture(por=0.30,sat=1.0,org=0.0), #(θwi=0.30,θm=0.70,θo=0.0,ϕ=0.30),
)
tempprofile_linear = TemperatureProfile(
    0.0u"m" => -30.0u"°C",
    10.0u"m" => -10.0u"°C", 
    1000.0u"m" => 10.2u"°C"
)
modelgrid = Grid(vcat(-2.0u"m":0.02u"m":-0.02u"m", CryoGrid.Presets.DefaultGrid_2cm))
z_top = -2.0u"m"
z_sub = map(knot -> knot.depth, soilprofile)
z_bot = modelgrid[end]
upperbc = TemperatureGradient(tair, NFactor())
initT = initializer(:T, tempprofile_linear)
@info "Building stratigraphy"
heatop = Heat.EnthalpyImplicit()
strat = @Stratigraphy(
    z_top => Top(upperbc),
    z_sub[1] => :lake => Lake(HeatBalance(heatop)),
    #z_sub[1] => :topsoil1 => Soil(HeatBalance(heatop), para=soilprofile[1].value),
    #z_sub[2] => :topsoil2 => Soil(HeatBalance(heatop), para=soilprofile[2].value),
    #z_sub[3] => :sediment1 => Soil(HeatBalance(heatop), para=soilprofile[3].value),
    #z_sub[4] => :sediment2 => Soil(HeatBalance(heatop), para=soilprofile[4].value),
    #z_sub[5] => :sediment3 => Soil(HeatBalance(heatop), para=soilprofile[5].value),
    z_bot => Bottom(GeothermalHeatFlux(0.053u"W/m^2"))
);
@info "Building tile"
tile = @time Tile(strat, modelgrid, initT)
# define time span, 5 years
tspan = (DateTime(2010,12,30), DateTime(2015,12,30))
tspan_sol = convert_tspan(tspan)
u0, du0 = @time initialcondition!(tile, tspan);
prob = CryoGridProblem(tile, u0, tspan, saveat=24*3600.0, savevars=(:θw,:T,))
lake_state = getstate(:lake, tile, u0, du0, prob.tspan[1])
# debug one step
CryoGrid.debug(true)
tile(du0, u0, prob.p, prob.tspan[1])
integrator = init(prob, LiteImplicitEuler(), dt=24*3600)
step!(integrator)
@info "Running model"
sol = @time solve(prob, LiteImplicitEuler(), dt=24*3600)
out = CryoGridOutput(sol)

# Plot the results
zs = [5,10,15,20,25,30,40,50,100,500,1000,5000]u"cm"
cg = Plots.cgrad(:copper,rev=true);
plot(out.T[Z(Near(zs))] |> ustrip, color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", title="", leg=false, dpi=150)
