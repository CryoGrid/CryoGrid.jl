using CryoGrid
const gridvals = vcat([0:0.02:2...,2.05:0.05:4.0...,
	4.1:0.1:10...,10.2:0.2:20...,21:1:30...,
	35:5:50...,60:10:100...,200:100:1000...]...)
# soil profile: depth => (total water, liquid water, mineral organic, porosity)
soilprofile = SoilProfile(
	0.0u"m" => (0.80,0.0,0.05,0.15,0.80),
	0.1u"m" => (0.80,0.0,0.15,0.05,0.80),
	0.4u"m" => (0.80,0.0,0.15,0.05,0.55),
	3.0u"m" => (0.50,0.0,0.50,0.0,0.50),
	10.0u"m" => (0.30,0.0,0.70,0.0,0.30),
)
tempprofile = TempProfile(
	0.0u"m" => 0.0u"°C",
	2.0u"m" => -1.0u"°C",
	5.0u"m" => -3.0u"°C",
	10.0u"m" => -6.0u"°C",
	25.0u"m" => -9.0u"°C",
	100.0u"m" => -9.0u"°C",
	1000.0u"m" => 10.2*u"°C"
)
strat = Stratigraphy(
	-2.0u"m" => Top(ConstantAirTemp(5.0u"°C")),
	0.0u"m" => Ground(:soil, Soil{Sand}(soilprofile), Heat{UT"J"}(tempprofile)),
	1000.0u"m" => Bottom(GeothermalHeatFlux(0.05u"J/s"))
)
grid = Grid(gridvals)
model = CryoGridSetup(strat,grid)
# define time span
tspan = [0.0u"s",180u"d"] |> ustrip
# CryoGrid front-end for ODEProblem
prob, out = CryoGridProblem(model,tspan,savevars=(soil=(:T,:θl),))
# solve discretized system w/ trapezoid method (Crank-Nicolson); save at 24 hour intervals
sol = @time solve(prob, Trapezoid(autodiff=false), abstol=0.1, saveat=24*3600.0)
# 6.516993 seconds (673.82 k allocations: 297.257 MiB, 0.65% gc time)

using Plots
Plots.pyplot()
Hgrid = round.(model.state.soil.grids.H,digits=2)
t = uconvert.(u"d",(sol.t)u"s")
zs = [1:10...,20:30:278...]
C(g::Plots.ColorGradient) = reshape(Plots.RGB[g[z] for z=LinRange(0,1,length(zs))],(1,:))
plot(ustrip(t),sol[zs,:]', label=Hgrid[zs]', color=Plots.cgrad(:berlin,rev=true) |> C,
	xlabel="Days", ylabel="H", leg=false, cbar=true)

T = hcat(map(x->x[1],out.saveval)...)
θl = hcat(map(x->x[2],out.saveval)...)
t = uconvert.(u"d",(out.t)u"s")
C(g::Plots.ColorGradient) = reshape(Plots.RGB[g[z] for z=LinRange(0,1,length(zs))],(1,:))
plot(ustrip(t), T[zs,:]'.-273.15, label=Hgrid[zs]', xlabel="Days", ylabel="Temperature (K)", c=Plots.cgrad(:berlin,rev=true) |> C, leg=false, cbar=true)
plot(ustrip(t), θl[zs,:]', label=Hgrid[zs]', xlabel="Days", ylabel="Liquid water content", leg=false, cbar=true)

using Interpolations
using MAT
interp = interpolate((Hgrid,out.t),T,Gridded(Linear()))
Tout = interp(Hgrid,sol.t)
plot(sol.t, Tout[zs,:]')
matwrite("/tmp/cryogrid_temps_1year.mat", Dict("data"=>Tout))
