Surface.surfaceproperties(seb::SurfaceEnergyBalance, sub::TestGroundLayer) = seb.para.soil

@testset "Surface Energy Balance" begin
    const_Tair = ConstantForcing(0.0u"°C", :Tair)
    const_pr = ConstantForcing(upreferred(1.0u"atm"), :pr)
    const_qh = ConstantForcing(0.001u"kg/kg", :qh)
    const_wind = ConstantForcing(1.0u"m/s", :wind)
    const_Lin = ConstantForcing(upreferred(300.0u"W/m^2"), :Lin)
    const_Sin = ConstantForcing(upreferred(200.0u"W/m^2"), :Sin)
    z = 2.0u"m"
    @testset "Iterative" begin
        seb = SurfaceEnergyBalance(
            const_Tair,
            const_pr,
            const_qh,
            const_wind,
            const_Lin,
            const_Sin,
            z,
            solscheme = Surface.Iterative(),
        )
        # TODO: would be good to have SEB work correctly with units to verify correctness;
        # however, at the moment it seems there are some constants that are used without units
        # which means that it unfortunately will not work
        seb = pstrip(seb)
        grid = Grid([0.0,0.1]u"m")
        toplayer = Top(seb)
        sublayer = TestGroundLayer(nothing)
        stop = Diagnostics.build_dummy_state(grid, toplayer, with_units=false)
        CryoGrid.initialcondition!(toplayer, stop)
        ssoil = (T=[-1.0],)
        sebstate = Surface.SEBState(seb, sublayer, stop, ssoil)
        seb_output = seb(sebstate)
        @test isfinite(seb_output.Qg)
        @test seb_output.Qg > zero(seb_output.Qg)
    end
    @testset "Analytical" begin
        seb = SurfaceEnergyBalance(
            const_Tair,
            const_pr,
            const_qh,
            const_wind,
            const_Lin,
            const_Sin,
            z,
            solscheme = Surface.Analytical(),
        )
        seb = pstrip(seb)
        grid = Grid([0.0,0.1]u"m")
        toplayer = Top(seb)
        sublayer = TestGroundLayer(nothing)
        stop = Diagnostics.build_dummy_state(grid, toplayer, with_units=false)
        CryoGrid.initialcondition!(toplayer, stop)
        ssoil = (T=[-1.0],)
        sebstate = Surface.SEBState(seb, sublayer, stop, ssoil)
        seb_output = seb(sebstate)
        @test isfinite(seb_output.Qg)
        @test seb_output.Qg > zero(seb_output.Qg)
    end
    @testset "Numerical" begin
        seb = SurfaceEnergyBalance(
        const_Tair,
        const_pr,
        const_qh,
        const_wind,
        const_Lin,
        const_Sin,
        z,
        solscheme=Surface.Numerical()
    )
        seb = pstrip(seb)
        grid = Grid([0.0,0.1]u"m")
        toplayer = Top(seb)
        sublayer = TestGroundLayer(nothing)
        stop = Diagnostics.build_dummy_state(grid, toplayer, with_units=false)
        CryoGrid.initialcondition!(toplayer, stop)
        ssoil = (T=[-1.0],)
        sebstate = Surface.SEBState(seb, sublayer, stop, ssoil)
        seb_output = solve(seb, sebstate)
        @test isfinite(seb_output.Qg)
        @test seb_output.Qg > zero(seb_output.Qg)
    end
end
