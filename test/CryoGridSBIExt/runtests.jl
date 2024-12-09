using CryoGrid
using OrdinaryDiffEq
using SimulationBasedInference
using Test

function cryogrid_test_problem(tspan)
    tempprofile = TemperatureProfile(
        0.0u"m" => -20.0u"°C",
        1000.0u"m" => 10.0u"°C",
    );
    
    # Here we use a simple single layer model with default soil parameters (50% porosity, no organic).
    soilprofile = SoilProfile(
        0.0u"m" => SimpleSoil()
    );
    heatop = Heat.Diffusion1D(:H)
    initT = initializer(:T, tempprofile)
    tile = CryoGrid.SoilHeatTile(
        heatop,
        ## 10 W/m^2 in and out, i.e. net zero flux
        PeriodicBC(HeatBalance, CryoGrid.Dirichlet, 365.0u"d", 15.0u"°C", 0.0, -5.0u"°C"),
        ConstantBC(HeatBalance, CryoGrid.Neumann, -0.1u"W/m^2"),
        soilprofile,
        inputs(),
        initT;
        grid=CryoGrid.DefaultGrid_10cm,
    );
    u0, _ = initialcondition!(tile, tspan);
    prob = CryoGridProblem(tile, u0, tspan)
    return prob
end

@testset "CryoGrid observables" begin
    tspan = (DateTime(2010,1,1), DateTime(2012,1,1))
    prob = cryogrid_test_problem(tspan)
    tile = Tile(prob.f)
    T_profile_observable = TemperatureProfileObservable(:Ts, [0.1,0.5,1.0,5.0,10.0], tspan, Day(1))
    T_var_observable = LayerVarObservable(:Ts, :ground1, :T, cells(tile.grid), tspan, Day(1))
    alt_profile_observable = ActiveLayerThicknessObservable(:alt, tspan)
    forward_prob = SimulatorForwardProblem(prob, T_profile_observable, alt_profile_observable)
    forward_sol = solve(forward_prob, Euler(), dt=120.0)
    model_Ts = SBI.getvalue(T_profile_observable)
    model_alt = SBI.getvalue(alt_profile_observable)
    @test size(model_Ts) == (5,365*2)
    @test size(model_alt) == (2,)
end
