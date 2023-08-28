@testset "CGEuler" begin
    @testset "water flow, bucket scheme" begin
        prob = test_water_flow_bucket_scheme()
        sol = solve(prob, CGEuler())
        out = CryoGridOutput(sol)
        # check that last grid cell is near saturation
        @test 1.0 - out.sat[end,end] < 0.001
    end
end
