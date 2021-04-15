using CryoGrid
using Test
using ForwardDiff
using ComponentArrays

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
            γ = 0.1
            f = McKenzie()
            sfcc = SFCC(f, SFCCNewtonSolver(tol=tol, onfail=:error))
            soil = Soil{Sand}(testprofile)
            heat = Heat{u"J"}(freezecurve=sfcc)
            L = heat.params.L
            @testset "Left tail" begin
                # set up single-grid-cell state vars
                θw,θl,θm,θo,θp = map(x -> [x], testprofile[1,:]) # convert to arrays
                T = [-5.0 + 273.15] # convert to K
                θl = f.(T,γ,θw,θp) # set liquid water content according to freeze curve
                C = heatcapacity.(soil.hcparams,θw,θl,θm,θo)
                H = enthalpy.(T.+1,C,L,θl) # compute enthalpy at +1 degree
                params = (γ=γ,)
                state = (T=T,C=C,H=H,θw=θw,θl=θl,θm=θm,θo=θo,θp=θp,params=params)
                sfcc(soil, heat, state)
                @test all(abs.((T.-273.15).-(H .- L.*θl)./C) .<= tol)
            end
            @testset "Right tail" begin
                # set up single-grid-cell state vars
                θw,θl,θm,θo,θp = map(x -> [x], testprofile[1,:]) # convert to arrays
                T = [5.0 + 273.15] # convert to K
                θl = f.(T,γ,θw,θp) # set liquid water content according to freeze curve
                C = heatcapacity.(soil.hcparams,θw,θl,θm,θo)
                H = enthalpy.(T.+1,C,L,θl) # compute enthalpy at +1 degree
                params = (γ=γ,)
                state = (T=T,C=C,H=H,θw=θw,θl=θl,θm=θm,θo=θo,θp=θp,params=params)
                sfcc(soil, heat, state)
                @test all(abs.((T.-273.15).-(H .- L.*θl)./C) .<= tol)
            end
            @testset "Near zero" begin
                # set up single-grid-cell state vars
                θw,θl,θm,θo,θp = map(x -> [x], testprofile[1,:]) # convert to arrays
                T = [-0.05 + 273.15] # convert to K
                θl = f.(T,γ,θw,θp) # set liquid water content according to freeze curve
                C = heatcapacity.(soil.hcparams,θw,θl,θm,θo)
                H = enthalpy.(T.+0.02,C,L,θl) # compute enthalpy at +.02 degree
                params = (γ=γ,)
                state = (T=T,C=C,H=H,θw=θw,θl=θl,θm=θm,θo=θo,θp=θp,params=params)
                sfcc(soil, heat, state)
                @test all(abs.((T.-273.15).-(H .- L.*θl)./C) .<= tol)
            end
        end
    end
    # TODO: DRY violation; a lot of this code is redundant and could possibly be
    # shared with the McKenzie freeze curve test.
    @testset "Van Genuchten freeze curve" begin
        @testset "Sanity checks" begin
            f = VanGenuchten()
            θsat = 0.8
            α = 4.0
            n = 2.0
            L = Heat{u"J"}().params.L
            @test isapprox(f(263.15,α,n,273.15,θsat,θsat,L), 0.0, atol=1e-6)
            @test f(273.15,α,n,273.15,θsat,θsat,L) ≈ θsat
            θl = f(273.05,α,n,273.15,θsat,θsat,L)
            @test θl > 0.0 && θl < 1.0
        end
        @testset "Newton solver checks" begin
            tol = 0.01
            θsat = 0.8
            α = 4.0
            n = 2.0
            Tₘ = 273.15
            f = VanGenuchten()
            sfcc = SFCC(f, SFCCNewtonSolver(tol=tol, onfail=:error))
            soil = Soil{Sand}(testprofile)
            heat = Heat{u"J"}(freezecurve=sfcc)
            L = heat.params.L
            @testset "Left tail" begin
                # set up single-grid-cell state vars
                θw,θl,θm,θo,θp = map(x -> [x], testprofile[1,:]) # convert to arrays
                T = [-5.0 + 273.15] # convert to K
                θl = f.(T,α,n,Tₘ,θw,θp,L) # set liquid water content according to freeze curve
                C = heatcapacity.(soil.hcparams,θw,θl,θm,θo)
                H = enthalpy.(T.+1,C,L,θl) # compute enthalpy at +1 degree
                params = (α=α,n=n,Tₘ=Tₘ)
                state = (T=T,C=C,H=H,θw=θw,θl=θl,θm=θm,θo=θo,θp=θp,params=params)
                @inferred sfcc(soil, heat, state)
                @test all(abs.((T.-273.15).-(H .- L.*θl)./C) .<= tol)
            end
            @testset "Right tail" begin
                # set up single-grid-cell state vars
                θw,θl,θm,θo,θp = map(x -> [x], testprofile[1,:]) # convert to arrays
                T = [5.0 + 273.15] # convert to K
                θl = f.(T,α,n,Tₘ,θw,θp,L) # set liquid water content according to freeze curve
                C = heatcapacity.(soil.hcparams,θw,θl,θm,θo)
                H = enthalpy.(T.+1,C,L,θl) # compute enthalpy at +1 degree
                params = (α=α,n=n,Tₘ=Tₘ)
                state = (T=T,C=C,H=H,θw=θw,θl=θl,θm=θm,θo=θo,θp=θp,params=params)
                @inferred sfcc(soil, heat, state)
                @test all(abs.((T.-273.15).-(H .- L.*θl)./C) .<= tol)
            end
            @testset "Near zero" begin
                # set up single-grid-cell state vars
                θw,θl,θm,θo,θp = map(x -> [x], testprofile[1,:]) # convert to arrays
                T = [-0.05 + 273.15] # convert to K
                θl = f.(T,α,n,Tₘ,θw,θp,L) # set liquid water content according to freeze curve
                C = heatcapacity.(soil.hcparams,θw,θl,θm,θo)
                H = enthalpy.(T.+0.02,C,L,θl) # compute enthalpy at +.02 degree
                params = (α=α,n=n,Tₘ=Tₘ)
                state = (T=T,C=C,H=H,θw=θw,θl=θl,θm=θm,θo=θo,θp=θp,params=params)
                @inferred sfcc(soil, heat, state)
                @test all(abs.((T.-273.15).-(H .- L.*θl)./C) .<= tol)
            end
        end
    end
    @testset "Newton solver autodiff" begin
        # set up
        tol = 0.01
        θsat = 0.8
        γ = 0.1
        f = McKenzie()
        sfcc = SFCC(f, SFCCNewtonSolver(tol=tol, onfail=:error))
        soil = Soil{Sand}(testprofile)
        heat = Heat{u"J"}(freezecurve=sfcc)
        L = heat.params.L
        θw,θl,θm,θo,θp = map(x -> [x], testprofile[1,:]) # convert to arrays
        T = [-0.1 + 273.15] # convert to K
        θl = f.(T,γ,θw,θp) # set liquid water content according to freeze curve
        C = heatcapacity.(soil.hcparams,θw,θl,θm,θo)
        H = enthalpy.(T.+0.09,C,L,θl) # compute enthalpy at +1 degree
        # test gradients
        p = ComponentArray(γ=γ)
        ∂f∂p = ForwardDiff.gradient(p ->  sum(f.(T,p.γ,θw,θp)), p)
        @test all(isfinite.(∂f∂p))
        function F(p)
            T = similar(T,eltype(p))
            T .= T
            C = similar(C,eltype(p))
            θl = similar(θl,eltype(p))
            state = (T=T,C=C,H=H,θw=θw,θl=θl,θm=θm,θo=θo,θp=θp,params=p)
            sfcc(soil, heat, state)
            state.T[1]
        end
        p = ComponentArray(γ=[γ])
        ∂F∂p = ForwardDiff.gradient(F, p)
        @test all(isfinite.(∂F∂p))
    end
end;

using BenchmarkTools
function benchmarksfcc()
    tol = 0.01
    θsat = 0.8
    α = 4.0
    n = 2.0
    Tₘ = 273.15
    f = VanGenuchten()
    sfcc = SFCC(f, SFCCNewtonSolver(tol=tol, α₀=1.0, τ=0.75, onfail=:error))
    soil = Soil{Sand}(testprofile)
    heat = Heat{u"J"}(freezecurve=sfcc)
    L = heat.params.L
    # set up multi-grid-cell state vars
    T = [-0.05 + 273.15 for i in 1:250] # convert to K
    θw,θl,θm,θo,θp = map(x -> x*ones(length(T)), testprofile[1,:]) # convert to arrays
    θl = f.(T,α,n,Tₘ,θw,θp,L) # set liquid water content according to freeze curve
    C = heatcapacity.(soil.hcparams,θw,θl,θm,θo)
    H = enthalpy.(T.+0.04,C,L,θl) # compute enthalpy at +.04 degree
    params = (α=α,n=n,Tₘ=Tₘ)
    state = (T=T,C=C,H=H,θw=θw,θl=θl,θm=θm,θo=θo,θp=θp,params=params)
    # sfcc.solver(soil, heat, state, sfcc.f, sfcc.∇f)
    # @time begin
    #     state.T .= -0.05 + 273.15
    #     sfcc.solver(soil, heat, state, sfcc.f, sfcc.∇f)
    # end
    @benchmark begin
        $state.T .= -0.05 + 273.15
        $sfcc($soil, $heat, $state)
    end
end
