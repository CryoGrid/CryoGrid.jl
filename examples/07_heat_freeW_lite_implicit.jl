# # [Fast heat conduction with CryoGridLite](@id example7)
# This example is very similar to [Example 1](@ref) but uses the
# fast implicit CryoGridLite solver of Langer et al. 2023.

using CryoGrid
using CryoGrid.LiteImplicit

forcings = loadforcings(CryoGrid.Presets.Forcings.Samoylov_ERA_MkL3_CCSM4_long_term);
tempprofile_linear = TemperatureProfile(
    0.0u"m" => -30.0u"°C",
    10.0u"m" => -10.0u"°C", 
    1000.0u"m" => 10.2u"°C"
)
z_top = -2.0u"m"
z_bot = 1000.0u"m"
upperbc = TemperatureGradient(forcings.Tair, NFactor())
initT = initializer(:T, tempprofile_linear)
@info "Building stratigraphy"
heatop = Heat.EnthalpyImplicit()
strat = @Stratigraphy(
    z_top => Top(upperbc),
    0.0u"m" => :topsoil1 => SimpleSoil(MineralOrganic(por=0.80,sat=1.0,org=0.75), heat=HeatBalance(heatop)),
    0.1u"m" => :topsoil2 => SimpleSoil(MineralOrganic(por=0.80,sat=1.0,org=0.25), heat=HeatBalance(heatop)),
    0.4u"m" => :sediment1 => SimpleSoil(MineralOrganic(por=0.55,sat=1.0,org=0.25), heat=HeatBalance(heatop)),
    3.0u"m" => :sediment2 => SimpleSoil(MineralOrganic(por=0.50,sat=1.0,org=0.0), heat=HeatBalance(heatop)),
    10.0u"m" => :sediment3 => SimpleSoil(MineralOrganic(por=0.30,sat=1.0,org=0.0), heat=HeatBalance(heatop)),
    z_bot => Bottom(GeothermalHeatFlux(0.053u"W/m^2"))
);
@info "Building tile"
modelgrid = CryoGrid.Presets.DefaultGrid_2cm
tile = Tile(strat, modelgrid, initT)
# Since the solver can take daily timesteps, we can easily specify longer simulation time spans at minimal cost.
# Here we specify a time span of 5 years.
tspan = (DateTime(2005,12,30), DateTime(2010,12,30))
tspan_sol = convert_tspan(tspan)
u0, du0 = initialcondition!(tile, tspan);
prob = CryoGridProblem(tile, u0, tspan, saveat=[], savevars=(:T,))
@info "Running model"
sol = @time solve(prob, LiteImplicitEuler())
out = CryoGridOutput(sol)

# Plot the results!
import Plots

zs = [5,10,15,20,25,30,40,50,100,500,1000,5000]u"cm"
cg = Plots.cgrad(:copper,rev=true);
Plots.plot(out.T[Z(Near(zs))], color=cg[LinRange(0.0,1.0,length(zs))]', ylabel="Temperature", title="", leg=false, dpi=150)
