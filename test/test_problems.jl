using CryoGrid

function test_heat_conduction_freeW_periodic_bc()
    # Select default grid and initial temperature profile.
    tempprofile = TemperatureProfile(
        0.0u"m" => -10.0u"°C",
        1000.0u"m" => 10.0u"°C",
    );
    soilprofile = SoilProfile(
        0.0u"m" => MineralOrganic(por=0.50),
    );
    heatop = Heat.EnthalpyForm(SFCCPreSolver())
    initT = initializer(:T, tempprofile)
    # Define the simulation time span.
    tspan = (0.0,2*24*3600.0)
    discretization = CryoGrid.Presets.DefaultGrid_10cm
    Qbot = -0.01u"W/m^2"
    tile = CryoGrid.Presets.SoilHeatTile(
        heatop,
        PeriodicBC(HeatBalance, CryoGrid.Neumann, 24*3600.0, 20.0, 0.0, -5.0),
        ConstantBC(HeatBalance, CryoGrid.Neumann, Qbot),
        soilprofile,
        initT;
        discretization,
    );
    u0, du0 = initialcondition!(tile, tspan);
    prob = CryoGridProblem(tile, u0, tspan, saveat=3600.0, savevars=(:T,:H))
    return prob
end

function test_water_flow_bucket_scheme(topflux=0.0u"m/s")
    initsat = initializer(:sat, 0.5)
    # Define the simulation time span.
    tspan = (0.0, 24*3600.0)
    discretization = CryoGrid.Presets.DefaultGrid_10cm
    strat = Stratigraphy(
        0.0u"m" => Top(ConstantBC(WaterBalance, CryoGrid.Neumann, topflux),),
        0.0u"m" => :ground => Ground(MineralOrganic(por=0.5, sat=0.5), heat=nothing, water=WaterBalance(BucketScheme())),
        1.0u"m" => Bottom(ConstantBC(WaterBalance, CryoGrid.Neumann, zero(topflux)),)
    )
    tile = Tile(strat, discretization, initsat)
    u0, du0 = initialcondition!(tile, tspan);
    prob = CryoGridProblem(tile, u0, tspan, saveat=60.0)
    return prob
end
