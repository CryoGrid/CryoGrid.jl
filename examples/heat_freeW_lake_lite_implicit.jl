using CryoGrid
using CryoGrid.LiteImplicit
using Plots

Plots.plotly()

CryoGrid.debug(true)

raw_forcings = loadforcings(CryoGrid.Forcings.Samoylov_ERA_MkL3_CCSM4_long_term);
Tair = raw_forcings.data.Tair
Ptot = uconvert.(u"m/s", raw_forcings.data.Ptot)
rainfall = Ptot.*(Tair .> 0u"°C")
snowfall = Ptot.*(Tair .<= 0u"°C")
forcings = rebuild(raw_forcings; Tair, rainfall, snowfall);
soilprofile = SoilProfile(
    0.0u"m" => SimpleSoil(por=0.80,sat=0.9,org=0.75),
    0.1u"m" => SimpleSoil(por=0.80,sat=1.0,org=0.25),
    0.4u"m" => SimpleSoil(por=0.55,sat=1.0,org=0.25),
    3.0u"m" => SimpleSoil(por=0.50,sat=1.0,org=0.0),
    10.0u"m" => SimpleSoil(por=0.30,sat=1.0,org=0.0),
)
tempprofile_linear = TemperatureProfile(
    -2.0u"m" => 0.0u"°C",
    0.0u"m" => 0.0u"°C",
    10.0u"m" => -10.0u"°C", 
    1000.0u"m" => 10.2u"°C"
)
modelgrid = Grid(vcat(-1.0u"m":0.02u"m":-0.02u"m", CryoGrid.DefaultGrid_2cm))
z_top = -1.0u"m"
z_sub = keys(soilprofile)
z_bot = modelgrid[end]
upperbc = TemperatureBC(Input(:Tair), NFactor(nf=0.5))
initT = initializer(:T, tempprofile_linear)
@info "Building stratigraphy"
heatop = Heat.EnthalpyImplicit()
strat = @Stratigraphy(
    z_top => Top(upperbc),
    -1.0u"m" => Lake(heat=HeatBalance(heatop)),
    0.0u"m" => Ground(soilprofile[1].value, heat=HeatBalance(heatop)),
    z_bot => Bottom(GeothermalHeatFlux(0.053u"W/m^2"))
);
@info "Building tile"
tile = @time Tile(strat, modelgrid, forcings, initT)
# define time span, 5 years
tspan = (DateTime(2010,12,30), DateTime(2012,12,30))
tspan_sol = convert_tspan(tspan)
u0, du0 = @time initialcondition!(tile, tspan);
prob = CryoGridProblem(tile, u0, tspan, saveat=24*3600.0, savevars=(:θw,:T,))
# set up integrator
integrator = init(prob, LiteImplicitEuler(), dt=24*3600)
# debug one step
step!(integrator)

for i in integrator
    state = getstate(integrator)
    if state.top.T_ub[1] > 0.0
        break
    end
end

step!(integrator)
state = getstate(integrator)
state.top.T_ub
state.lake.T

@info "Running model"
sol = @time solve(prob, LiteImplicitEuler(), dt=24*3600.0)
out = CryoGridOutput(sol)

# Plot the results
zs = [-100,-50.0,-10.0,11.0]u"cm"
cg = Plots.cgrad(:copper,rev=true);
plot(convert_t.(dims(out.T, Ti)), Array(out.T[Z(Near(zs))])' |> ustrip, color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", title="", dpi=150)
plot!(convert_t.(dims(out.T, Ti)), t -> forcings.Tair.(t), c=:blue, linestyle=:dash)
