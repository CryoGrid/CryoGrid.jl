using CryoGrid
using CryoGrid.LiteImplicit
using Dates
using Plots
Plots.default(tickfont=(10,"Dejavu Sans"), titlefont=(12,"Dejavu Sans"))

function heat_conduction_linear_periodic_ub(T₀, A, P, α)
    T(z,t) = T₀ + A*exp(-z*sqrt(π/(α*P)))*sin(2π*t/P - z*sqrt(π/(α*P)))
    return T
end

P = 365*24*3600.0 # 1 year
A = 1.0 # 1°C amplitude
T₀ = 1.0 # average temperature of -1.0°C

# @testset "CryoGridLite" begin
z_top = 0.0u"m"
z_bot = 1000.0u"m"
heatop = Heat.EnthalpyImplicit()
# heatop = Heat.InverseEnthalpy()
soil = Soil(HomogeneousMixture(por=0.0, org=0.0), heat=HeatBalance(heatop))
strat = @Stratigraphy(
    z_top => Top(PeriodicBC(HeatBalance, CryoGrid.Dirichlet, P, 1.0, 0.0, T₀)),
    z_top => :soil => soil,
    z_bot => Bottom(ConstantFlux(HeatBalance, 0.0))
);
α = upreferred(strat.soil.prop.heat.kh_m) / upreferred(strat.soil.prop.heat.ch_m)
T_analytic = heat_conduction_linear_periodic_ub(T₀, A, P, ustrip(α))
initT = initializer(:T, (layer, proc, state) -> state.T .= T_analytic.(cells(state.grid), 0.0))
modelgrid = CryoGrid.Presets.DefaultGrid_2cm
# modelgrid = Grid(z_top:0.02u"m":z_bot)
# modelgrid = Grid(vcat(0.0:0.02:1.0, 1.05:0.05:5.0, 5.1:0.1:10.0)*u"m")
tile = Tile(strat, modelgrid, initT)
# define time span, 5 years
tspan = (0.0,5*365*24*3600.0)
u0, du0 = initialcondition!(tile, tspan);
T0 = getvar(:T, tile, u0; interp=false)
prob = CryoGridProblem(tile, u0, tspan, saveat=24*3600.0, savevars=(:T,:kc,:C))
sol = @time solve(prob, LiteImplicitEuler(), dt=24*3600)
# sol = @time solve(prob, ImplicitEuler())
out = CryoGridOutput(sol)
Ts = out.T[Z(0.0u"m"..10.0u"m")]
ts = ustrip.(u"d", (tspan[1]:24*3600:tspan[end])*u"s")
zs = collect(dims(Ts, Z))
p1 = heatmap(ts, ustrip.(zs), (t,z) -> ustrip(Ts[Z(At(z*u"m")), Ti(At(convert_t(t*24*3600)))]), yflip=true, title="CryoGridLite", xlabel="Days", ylabel="Depth (m)", cbar=false)
# check analytical solution
p2 = heatmap(ts, zs, (t,z) -> T_analytic(ustrip(z),t*24*3600), yflip=true, title="Analytical", cbar=true)
plot(p1, p2, size=(1200,500), dpi=300, layout=@layout([a{0.45w} b{0.55w}]))

plot(ts, t -> T_analytic(1.01, t*24*3600), leg=nothing)
plot!(ts, t -> ustrip(Ts[Z(Near(1.01u"m")),Ti(At(convert_t(t*24*3600)))]))

# Stefan
z_top = 0.0u"m"
z_bot = 1000.0u"m"
heatop = Heat.EnthalpyImplicit()
# heatop = Heat.InverseEnthalpy(SFCCPreSolver())
soil = Soil(HomogeneousMixture(por=0.3, sat=1.0, org=0.0), heat=HeatBalance(heatop))
strat = @Stratigraphy(
    z_top => Top(ConstantTemperature(1.0u"°C")),
    z_top => :soil => soil,
    z_bot => Bottom(ConstantFlux(HeatBalance, 0.0))
);
α = upreferred(strat.soil.prop.heat.kh_m) / upreferred(strat.soil.prop.heat.ch_m)
T_analytic = heat_conduction_linear_periodic_ub(T₀, A, P, ustrip(α))
initT = initializer(:T, -1.0)
modelgrid = CryoGrid.Presets.DefaultGrid_2cm
tile = Tile(strat, modelgrid, initT)
# define time span, 5 years
tspan = (0.0, 5*365*24*3600.0)
u0, du0 = initialcondition!(tile, tspan);
T0 = getvar(:T, tile, u0)
prob = CryoGridProblem(tile, u0, tspan, saveat=24*3600.0, savevars=(:T,:θw,:kc,:C))
sol = @time solve(prob, LiteImplicitEuler(), dt=24*3600)
# sol = @time solve(prob, SSPRK43(), reltol=1e-6, saveat=24*3600.0, progress=true)
out = CryoGridOutput(sol)
plot(out.T[Z(1:10)])

kh_w, kh_i, kh_a, kh_m, kh_o = Heat.thermalconductivities(strat.soil)
θ_s = (θw=0.0, θi=porosity(strat.soil), θa=0.0, θm=mineral(strat.soil), θo=organic(strat.soil))
θ_l = (θw=porosity(strat.soil), θi=0.0, θa=0.0, θm=mineral(strat.soil), θo=organic(strat.soil))
k_s = Heat.thermalconductivity(strat.soil, strat.soil.heat, θ_s...)
k_l = Heat.thermalconductivity(strat.soil, strat.soil.heat, θ_l...)
c_s = Heat.heatcapacity(strat.soil, strat.soil.heat, θ_s...)
c_l = Heat.heatcapacity(strat.soil, strat.soil.heat, θ_l...)
stefan_prob = StefanProblem(p=StefanParameters(T_s=-1.0u"°C", T_l=1.0u"°C"; k_s, k_l, c_s, c_l, θwi=0.3))
stefan_sol = solve(stefan_prob)
# heatmap(0.0:3600.0:12*30*24*3600.0, 0.0:0.001:0.5, (t,z) -> ustrip(stefan_sol(z*u"m", t*u"s")), yflip=true)
ts = ustrip.(u"d", (tspan[1]:24*3600:tspan[end])*u"s")
plot(ustrip.(ts), t -> ustrip(stefan_sol(uconvert(u"s",t*u"d"))), linewidth=2, ylabel="Thaw depth (m)", xlabel="Days", yflip=true, leg=:topright, c=:black, label="Analytical", dpi=300)
td = Array(Diagnostics.thawdepth(out.θw./0.3, modelgrid))#[1:end-1]
plot!(ustrip.(ts), ustrip.(td), label="Model", linewidth=2, linestyle=:dash)
savefig("output/cglite_stefan_comparison_soil.svg")

# end
