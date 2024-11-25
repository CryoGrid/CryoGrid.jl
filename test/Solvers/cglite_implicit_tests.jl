using CryoGrid
using CryoGrid.LiteImplicit
using Dates
using Test

function thawdepth(θw::AbstractDimArray, grid::Grid; Tmelt=0.0u"°C")
    θw = permutedims(θw, (Ti,Z))
    thick = Δ(edges(grid))
    return DimArray(sum(θw.*reshape(thick, 1, length(thick)), dims=Z)[:,1], (dims(θw,Ti),))
end

@testset "CryoGridLite" begin
    # Linear heat conduction
    @testset "Linear heat with periodic boundary" begin
        P = 365*24*3600.0 # 1 year
        A = 1.0 # 1°C amplitude
        T₀ = 1.0 # average temperature of -1.0°C
        z_top = 0.0u"m"
        z_bot = 1000.0u"m"
        heatop = Heat.EnthalpyImplicit()
        soil = Ground(SimpleSoil(por=0.0, org=0.0), heat=HeatBalance(heatop))
        strat = @Stratigraphy(
            z_top => Top(PeriodicBC(HeatBalance, CryoGrid.Dirichlet, P, 1.0, 0.0, T₀)),
            z_top => :soil => soil,
            z_bot => Bottom(ConstantFlux(HeatBalance, 0.0))
        );
        α = upreferred(strat.soil.para.heat.kh_m) / upreferred(strat.soil.para.heat.ch_m)
        T_analytic = Heat.heat_conduction_linear_periodic_ub(T₀, A, P, ustrip(α))
        initT = initializer(:T, (layer, state) -> state.T .= T_analytic.(cells(state.grid), 0.0))
        modelgrid = CryoGrid.DefaultGrid_2cm
        # modelgrid = Grid(z_top:0.02u"m":z_bot)
        # modelgrid = Grid(vcat(0.0:0.02:1.0, 1.05:0.05:5.0, 5.1:0.1:10.0)*u"m")
        tile = Tile(strat, modelgrid, initT)
        # define time span, 5 years
        tspan = (0.0,5*365*24*3600.0)
        u0, du0 = initialcondition!(tile, tspan);
        T0 = getvar(:T, tile, u0; interp=false)
        prob = CryoGridProblem(tile, u0, tspan, saveat=24*3600.0, savevars=(:T,:kc,:C))
        sol = solve(prob, LiteImplicitEuler(), dt=24*3600)
        out = CryoGridOutput(sol)
        Ts = out.T[Z(0.0u"m"..10.0u"m")]
        ts = ustrip.(u"d", (tspan[1]:24*3600:tspan[end])*u"s")
        zs = collect(dims(Ts, Z))
        err = [ustrip(Ts[Z(At(z*u"m")), Ti(At(convert_t(t*24*3600)))]) - T_analytic(ustrip(z),t*24*3600) for t in ts, z in ustrip.(zs)]
        @test all(abs.(err) .< 0.01)
    end

    # Stefan
    @testset "Two-phase Stefan solution" begin
        z_top = 0.0u"m"
        z_bot = 1000.0u"m"
        heatop = Heat.EnthalpyImplicit()
        # heatop = Heat.Diffusion1D(:H)
        soil = Ground(SimpleSoil(por=0.3, sat=1.0, org=0.0), heat=HeatBalance(heatop), water=nothing)
        strat = @Stratigraphy(
            z_top => Top(ConstantTemperature(1.0u"°C")),
            z_top => :soil => soil,
            z_bot => Bottom(ConstantFlux(HeatBalance, 0.0))
        );
        initT = initializer(:T, -1.0)
        modelgrid = CryoGrid.DefaultGrid_2cm
        tile = Tile(strat, modelgrid, initT)
        # define time span, 5 years
        tspan = (0.0, 5*365*24*3600.0)
        u0, du0 = initialcondition!(tile, tspan);
        T0 = getvar(:T, tile, u0)
        prob = CryoGridProblem(tile, u0, tspan, saveat=24*3600.0, savevars=(:T,:θw,:kc,:C))
        sol = solve(prob, LiteImplicitEuler(), dt=24*3600)
        out = CryoGridOutput(sol)

        θm = mineral(strat.soil)
        θo = organic(strat.soil)
        θ_s = (θw=0.0, θi=porosity(strat.soil), θa=0.0, θm=θm, θo=θo)
        θ_l = (θw=porosity(strat.soil), θi=0.0, θa=0.0, θm=θm, θo=θo)
        k_vals = Heat.thermalconductivities(strat.soil)
        c_vals = Heat.heatcapacities(strat.soil)
        thermalprops = Heat.thermalproperties(strat.soil)
        k_s = thermalprops.thermcond(k_vals, θ_s)
        k_l = thermalprops.thermcond(k_vals, θ_l)
        c_s = thermalprops.heatcap(c_vals, θ_s)
        c_l = thermalprops.heatcap(c_vals, θ_l)
        stefan_prob = StefanProblem(p=StefanParameters(T_s=-1.0u"°C", T_l=1.0u"°C"; k_s, k_l, c_s, c_l, θwi=0.3))
        stefan_sol = solve(stefan_prob)
        ts = ustrip.(u"d", (tspan[1]:24*3600:tspan[end])*u"s")
        td = Array(thawdepth(out.θw./0.3, modelgrid))
        td_true = stefan_sol.(uconvert.(u"s", ts.*u"d"))
        @test all(abs.(td .- td_true) .< 0.01u"m")
    end
end
