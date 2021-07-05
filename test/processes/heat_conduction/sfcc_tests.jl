using CryoGrid
using Test
using ForwardDiff
using ComponentArrays

testprofile = SoilProfile(
	0.0u"m" => SoilProperties(χ=0.0,ϕ=0.80,θ=1.0,ω=0.5),
	1.0u"m" => SoilProperties(χ=0.0,ϕ=0.80,θ=1.0,ω=0.5),
)
soilcomps = begin 
    comps = [CryoGrid.Layers.soilcomp(Val{var}(),testprofile[var=:χ],testprofile[var=:ϕ],testprofile[var=:θ],testprofile[var=:ω]) for var in [:θx,:θp,:θm,:θo]]
    reduce(hcat, [comps[1] .+ comps[2], comps[3], comps[4], comps[2]])
end

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
            soil = Soil(testprofile)
            heat = Heat{:H}(freezecurve=sfcc)
            L = heat.params.L
            @testset "Left tail" begin
                # set up single-grid-cell state vars
                θw,θm,θo,θp = map(x -> [x], soilcomps[1,:]) # convert to arrays
                T = [-5.0 + 273.15] # convert to K
                θl = f.(T,γ,θw,θp) # set liquid water content according to freeze curve
                C = heatcapacity.(soil.params,θw,θl,θm,θo)
                H = let T = T.+1.0,
                    θl = f.(T,γ,θw,θp)
                    C = heatcapacity.(soil.params,θw,θl,θm,θo);
                   enthalpy.(T,C,L,θl) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,Ceff=similar(C),H=H,θw=θw,θl=θl,θm=θm,θo=θo,θp=θp,params=(γ=γ,))
                sfcc(soil, heat, state)
                @test all(abs.((T.-273.15).-(H .- L.*θl)./C) .<= tol)
            end
            @testset "Right tail" begin
                # set up single-grid-cell state vars
                θw,θm,θo,θp = map(x -> [x], soilcomps[1,:]) # convert to arrays
                T = [5.0 + 273.15] # convert to K
                θl = f.(T,γ,θw,θp) # set liquid water content according to freeze curve
                C = heatcapacity.(soil.params,θw,θl,θm,θo)
                H = let T = T.-1.0,
                    θl = f.(T,γ,θw,θp)
                    C = heatcapacity.(soil.params,θw,θl,θm,θo);
                   enthalpy.(T,C,L,θl) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,Ceff=similar(C),H=H,θw=θw,θl=θl,θm=θm,θo=θo,θp=θp,params=(γ=γ,))
                sfcc(soil, heat, state)
                @test all(abs.((T.-273.15).-(H .- L.*θl)./C) .<= tol)
            end
            @testset "Near zero" begin
                # set up single-grid-cell state vars
                θw,θm,θo,θp = map(x -> [x], soilcomps[1,:]) # convert to arrays
                T = [-0.05 + 273.15] # convert to K
                θl = f.(T,γ,θw,θp) # set liquid water content according to freeze curve
                C = heatcapacity.(soil.params,θw,θl,θm,θo)
                H = let T = T.+0.04,
                    θl = f.(T,γ,θw,θp)
                    C = heatcapacity.(soil.params,θw,θl,θm,θo);
                   enthalpy.(T,C,L,θl) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,Ceff=similar(C),H=H,θw=θw,θl=θl,θm=θm,θo=θo,θp=θp,params=(γ=γ,))
                sfcc(soil, heat, state)
                @test all(abs.((T.-273.15).-(H .- L.*θl)./C) .<= tol)
            end
        end
    end
    # TODO: DRY violation; a lot of this code is redundant and could possibly be
    # shared between different freeze curve tests.
    @testset "Westermann freeze curve" begin
        @testset "Sanity checks" begin
            f = Westermann()
            θtot = 0.8
            δ = 0.1
            @test isapprox(f(263.15,δ,θtot), 0.0, atol=1e-2)
            @test f(273.15,δ,θtot) ≈ θtot
            θl = f(273.05,δ,θtot)
            @test θl > 0.0 && θl < 1.0
        end
        @testset "Newton solver checks" begin
            tol = 0.01
            δ = 0.1
            f = Westermann()
            sfcc = SFCC(f, SFCCNewtonSolver(tol=tol, onfail=:error))
            soil = Soil(testprofile)
            heat = Heat{:H}(freezecurve=sfcc)
            L = heat.params.L
            @testset "Left tail" begin
                # set up single-grid-cell state vars
                θw,θm,θo,θp = map(x -> [x], soilcomps[1,:]) # convert to arrays
                T = [-5.0 + 273.15] # convert to K
                θl = f.(T,δ,θw) # set liquid water content according to freeze curve
                C = heatcapacity.(soil.params,θw,θl,θm,θo)
                H = let T = T.+1.0,
                    θl = f.(T,δ,θw)
                    C = heatcapacity.(soil.params,θw,θl,θm,θo);
                   enthalpy.(T,C,L,θl) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,Ceff=similar(C),H=H,θw=θw,θl=θl,θm=θm,θo=θo,θp=θp,params=(δ=δ,))
                sfcc(soil, heat, state)
                @test all(abs.((T.-273.15).-(H .- L.*θl)./C) .<= tol)
            end
            @testset "Right tail" begin
                # set up single-grid-cell state vars
                θw,θm,θo,θp = map(x -> [x], soilcomps[1,:]) # convert to arrays
                T = [5.0 + 273.15] # convert to K
                θl = f.(T,δ,θw) # set liquid water content according to freeze curve
                C = heatcapacity.(soil.params,θw,θl,θm,θo)
                H = let T = T.-1,
                    θl = f.(T,δ,θw)
                    C = heatcapacity.(soil.params,θw,θl,θm,θo);
                   enthalpy.(T,C,L,θl) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,Ceff=similar(C),H=H,θw=θw,θl=θl,θm=θm,θo=θo,θp=θp,params=(δ=δ,))
                sfcc(soil, heat, state)
                @test all(abs.((T.-273.15).-(H .- L.*θl)./C) .<= tol)
            end
            @testset "Near zero" begin
                # set up single-grid-cell state vars
                θw,θm,θo,θp = map(x -> [x], soilcomps[1,:]) # convert to arrays
                T = [-0.05 + 273.15] # convert to K
                θl = f.(T,δ,θw) # set liquid water content according to freeze curve
                C = heatcapacity.(soil.params,θw,θl,θm,θo)
                H = let T = T.+0.04,
                    θl = f.(T,δ,θw)
                    C = heatcapacity.(soil.params,θw,θl,θm,θo);
                   enthalpy.(T,C,L,θl) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,Ceff=similar(C),H=H,θw=θw,θl=θl,θm=θm,θo=θo,θp=θp,params=(δ=δ,))
                sfcc(soil, heat, state)
                @test all(abs.((T.-273.15).-(H .- L.*θl)./C) .<= tol)
            end
        end
    end
    @testset "Van Genuchten freeze curve" begin
        @testset "Sanity checks" begin
            f = VanGenuchten()
            θsat = 0.8
            α = 4.0
            n = 2.0
            L = Heat{:H}().params.L
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
            soil = Soil(testprofile)
            heat = Heat{:H}(freezecurve=sfcc)
            L = heat.params.L
            @testset "Left tail" begin
                # set up single-grid-cell state vars
                θw,θm,θo,θp = map(x -> [x], soilcomps[1,:]) # convert to arrays
                T = [-5.0 + 273.15] # convert to K
                θl = f.(T,α,n,Tₘ,θw,θp,L) # set liquid water content according to freeze curve
                C = heatcapacity.(soil.params,θw,θl,θm,θo)
                H = let T = T.+1.0,
                    θl = f.(T,α,n,Tₘ,θw,θp,L)
                    C = heatcapacity.(soil.params,θw,θl,θm,θo);
                   enthalpy.(T,C,L,θl) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,Ceff=similar(C),H=H,θw=θw,θl=θl,θm=θm,θo=θo,θp=θp,params=(α=α,n=n,Tₘ=Tₘ))
                @inferred sfcc(soil, heat, state)
                @test all(abs.((T.-273.15).-(H .- L.*θl)./C) .<= tol)
            end
            @testset "Right tail" begin
                # set up single-grid-cell state vars
                θw,θm,θo,θp = map(x -> [x], soilcomps[1,:]) # convert to arrays
                T = [5.0 + 273.15] # convert to K
                θl = f.(T,α,n,Tₘ,θw,θp,L) # set liquid water content according to freeze curve
                C = heatcapacity.(soil.params,θw,θl,θm,θo)
                H = let T = T.-1.0,
                    θl = f.(T,α,n,Tₘ,θw,θp,L)
                    C = heatcapacity.(soil.params,θw,θl,θm,θo);
                   enthalpy.(T,C,L,θl) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,Ceff=similar(C),H=H,θw=θw,θl=θl,θm=θm,θo=θo,θp=θp,params=(α=α,n=n,Tₘ=Tₘ))
                @inferred sfcc(soil, heat, state)
                @test all(abs.((T.-273.15).-(H .- L.*θl)./C) .<= tol)
            end
            @testset "Near zero" begin
                # set up single-grid-cell state vars
                θw,θm,θo,θp = map(x -> [x], soilcomps[1,:]) # convert to arrays
                T = [-0.05 + 273.15] # convert to K
                θl = f.(T,α,n,Tₘ,θw,θp,L) # set liquid water content according to freeze curve
                C = heatcapacity.(soil.params,θw,θl,θm,θo)
                H = let T = T.+0.04,
                    θl = f.(T,α,n,Tₘ,θw,θp,L)
                    C = heatcapacity.(soil.params,θw,θl,θm,θo);
                   enthalpy.(T,C,L,θl) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,Ceff=similar(C),H=H,θw=θw,θl=θl,θm=θm,θo=θo,θp=θp,params=(α=α,n=n,Tₘ=Tₘ))
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
        soil = Soil(testprofile)
        heat = Heat{:H}(freezecurve=sfcc)
        L = heat.params.L
        θw,θm,θo,θp = map(x -> [x], soilcomps[1,:]) # convert to arrays
        T = [-0.1 + 273.15] # convert to K
        θl = f.(T,γ,θw,θp) # set liquid water content according to freeze curve
        C = heatcapacity.(soil.params,θw,θl,θm,θo)
        H = enthalpy.(T.+0.09,C,L,θl) # compute enthalpy at +1 degree
        # test gradients
        p = ComponentArray(γ=γ)
        ∂f∂p = ForwardDiff.gradient(p ->  sum(f.(T,p.γ,θw,θp)), p)
        @test all(isfinite.(∂f∂p))
        function F(p)
            T_ = similar(T,eltype(p))
            T_ .= T
            C = similar(C,eltype(p))
            θl = similar(θl,eltype(p))
            state = (T=T_,C=C,Ceff=similar(C),H=H,θw=θw,θl=θl,θm=θm,θo=θo,θp=θp,params=(γ=p.γ,))
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
    soil = Soil(testprofile)
    heat = Heat{:H}(freezecurve=sfcc)
    L = heat.params.L
    # set up multi-grid-cell state vars
    T = [-15.0 + 273.15 for i in 1:10] # convert to K
    θw,θl,θm,θo,θp = map(x -> x*ones(length(T)), testprofile[1,:]) # convert to arrays
    θl = f.(T,α,n,Tₘ,θw,θp,L) # set liquid water content according to freeze curve
    C = heatcapacity.(soil.params,θw,θl,θm,θo)
    H = let T = T.+14.999,
            θl = f.(T,α,n,Tₘ,θw,θp,L)
            C = heatcapacity.(soil.params,θw,θl,θm,θo);
        enthalpy.(T,C,L,θl) # compute enthalpy at "true" temperature
    end
    state = (T=T,C=C,Ceff=similar(C),H=H,θw=θw,θl=θl,θm=θm,θo=θo,θp=θp,params=(α=α,n=n,Tₘ=Tₘ))
    # sfcc.solver(soil, heat, state, sfcc.f, sfcc.∇f)
    # @time begin
    #     state.T .= -0.05 + 273.15
    #     sfcc.solver(soil, heat, state, sfcc.f, sfcc.∇f)
    # end
    result = @benchmark begin
        $sfcc($soil, $heat, $state)
    end
    show(stdout, "text/plain", result)
    println("\nsolution: $(state.T[1])")
end
