using CryoGrid
const gridvals = vcat([-2.0,0:0.02:2...,2.05:0.05:4.0...,
	4.1:0.1:10...,10.2:0.2:20...,21:1:30...,
	35:5:50...,60:10:100...,200:100:1000...]...)
# soil profile: depth => (total water, liquid water, mineral organic, porosity)
soilprofile = SoilProfile(
	-2.0u"m" => (0.0,0.0,0.0,0.0,0.0),
	0.0u"m" => (0.80,0.0,0.05,0.15,0.80),
	0.1u"m" => (0.80,0.0,0.15,0.05,0.80),
	0.4u"m" => (0.80,0.0,0.15,0.05,0.55),
	3.0u"m" => (0.50,0.0,0.50,0.0,0.50),
	10.0u"m" => (0.30,0.0,0.70,0.0,0.30),
)
tempprofile = TempProfile(
	-2.0u"m" => 5.0u"°C",
	0.0u"m" => 0.0u"°C",
	2.0u"m" => -2.0u"°C",
	5.0u"m" => -7.0u"°C",
	10.0u"m" => -9.0u"°C",
	25.0u"m" => -9.0u"°C",
	100.0u"m" => -8.0u"°C",
	1000.0u"m" => 10.2u"°C"
)
strat = Stratigraphy(
	-2.0U"m" => Top(ConstantAirTemp(5.0U"°C")),
	0.0U"m" => Ground(:soil, Soil{Sand}(soilprofile), Heat{UT"J"}(tempprofile)),
	1000.0U"m" => Bottom(GeothermalHeatFlux(0.05U"W/m^2"))
)
grid = Grid(gridvals)
model = CryoGridSetup(strat,grid)
u0, du0 = initialcondition!(model)
# define time span
tspan = [0.0u"s",365u"d"] |> Tuple
# setup ODEProblem
prob = ODEProblem(model,u0,ustrip.(tspan))
# solve discretized system; save at 3 hour intervals
sol = @time solve(prob,alg_hints=[:stiff], saveat=3.0*3600.0)

using Plots
Hgrid = round.(model.state.soil.grids.H,digits=2)
t = uconvert.(u"d",(sol.t)u"s")
zs = [1,5,10,30,80,120,150,200,250,270,278]
plot(ustrip(t),sol[zs,:]', label=Hgrid[zs]', xlabel="Days", ylabel="H",
	title="CryoGrid heat conduction, constant air temp.")

function testalloc(u0,du0)
	du0 .= 0.0
	model(du0,u0,nothing,0.0)
end

using BenchmarkTools
@benchmark testalloc($u0,$du0)
