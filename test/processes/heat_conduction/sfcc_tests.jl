using CryoGrid
using Test

testprofile = SoilProfile(
	0.0u"m" => (0.80,0.0,0.05,0.15,0.80),
	1.0u"m" => (0.80,0.0,0.05,0.15,0.80),
)

@testset "SFCC" begin
    @testset "McKenzie freeze curve" begin
        @testset "Sanity checks" begin
            f = McKenzie()
            θsat = 0.8
            @test isapprox(f(263.15,1,θsat,θsat), 0.0, atol=1e-6)
            @test f(273.15,1,θsat,θsat) ≈ θsat
            θl = f(273.05,1,θsat,θsat)
            @test θl > 0.0 && θl < 1.0
        end
        @testset "Newton solver checks" begin
            tol = 0.01
            δ = 0.1
            f = McKenzie()
            sfcc = SFCC(f, SFCCNewtonSolver(tol=tol, onfail=:error))
            soil = Soil{Sand}(testprofile)
            heat = Heat{u"J"}(freezecurve=sfcc)
            L = heat.params.L
            @testset "Left tail" begin
                # set up single-grid-cell state vars
                θw,θl,θm,θo,θp = map(x -> [x], testprofile) # convert to arrays
                T = [-5.0 + 273.15] # convert to K
                θl = f.(T,δ,θw,θp) # set liquid water content according to freeze curve
                C = heatcapacity.(soil.hcparams,θw,θl,θm,θo)
                H = enthalpy.(T.+1,C,L,θl) # compute enthalpy at +1 degree
                params = (δ=δ,)
                state = (T=T,C=C,H=H,θw=θw,θl=θl,θm=θm,θo=θo,θp=θp,params=params)
                sfcc(soil, heat, state)
                @test all(abs.((T.-273.15).-(H .- L.*θl)./C) .<= tol)
            end
            @testset "Right tail" begin
                # set up single-grid-cell state vars
                θw,θl,θm,θo,θp = map(x -> [x], testprofile) # convert to arrays
                T = [5.0 + 273.15] # convert to K
                θl = f.(T,δ,θw,θp) # set liquid water content according to freeze curve
                C = heatcapacity.(soil.hcparams,θw,θl,θm,θo)
                H = enthalpy.(T.+1,C,L,θl) # compute enthalpy at +1 degree
                params = (δ=δ,)
                state = (T=T,C=C,H=H,θw=θw,θl=θl,θm=θm,θo=θo,θp=θp,params=params)
                sfcc(soil, heat, state)
                @test all(abs.((T.-273.15).-(H .- L.*θl)./C) .<= tol)
            end
            @testset "Near zero" begin
                # set up single-grid-cell state vars
                θw,θl,θm,θo,θp = map(x -> [x], testprofile) # convert to arrays
                T = [-0.05 + 273.15] # convert to K
                θl = f.(T,δ,θw,θp) # set liquid water content according to freeze curve
                C = heatcapacity.(soil.hcparams,θw,θl,θm,θo)
                H = enthalpy.(T.+0.02,C,L,θl) # compute enthalpy at +.02 degree
                params = (δ=δ,)
                state = (T=T,C=C,H=H,θw=θw,θl=θl,θm=θm,θo=θo,θp=θp,params=params)
                sfcc(soil, heat, state)
                @test all(abs.((T.-273.15).-(H .- L.*θl)./C) .<= tol)
            end
        end
    end
end
