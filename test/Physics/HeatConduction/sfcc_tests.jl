using CryoGrid
using ComponentArrays
using ForwardDiff
using Setfield
using Test

@testset "SFCC" begin
    Tₘ = 0.0
    θres = 0.0
    soil = Soil(para=CharacteristicFractions()) |> stripparams |> stripunits
    θw = Soils.soilcomponent(Val{:θw}(), soil.para)
    θp = Soils.soilcomponent(Val{:θp}(), soil.para)
    θm = Soils.soilcomponent(Val{:θm}(), soil.para)
    θo = Soils.soilcomponent(Val{:θo}(), soil.para)
    @testset "McKenzie freeze curve" begin
        @testset "Sanity checks" begin
            f = McKenzie() |> stripparams |> stripunits
            θsat = 0.8
            @test isapprox(f(-10.0,Tₘ,θres,θsat,θsat,1.0), 0.0, atol=1e-6)
            @test f(0.0,Tₘ,θres,θsat,θsat,1) ≈ θsat
            θl = f(-0.1,Tₘ,θres,θsat,θsat,1)
            @test θl > 0.0 && θl < 1.0
        end
        @testset "Newton solver checks" begin
            abstol = 1e-2
            reltol = 1e-4
            γ = 0.1
            f = McKenzie() |> stripparams |> stripunits
            sfcc = SFCC(f, SFCCNewtonSolver(abstol=abstol, reltol=reltol, onfail=:error))
            heat = Heat(freezecurve=sfcc) |> stripparams |> stripunits
            L = heat.L
            @testset "Left tail" begin
                T = [-5.0]
                θl = f.(T,Tₘ,θres,θp,θw,γ) # set liquid water content according to freeze curve
                C = heatcapacity.(soil,heat,θw,θl,θm,θo)
                H = let T = T.+1.0,
                    θl = f.(T,Tₘ,θres,θp,θw,γ)
                    C = heatcapacity.(soil,heat,θw,θl,θm,θo);
                   enthalpy.(T,C,L,θl) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,dHdT=similar(C),H=H,θl=θl,)
                sfcc(soil, heat, state)
                @test all(abs.(T.-(H .- L.*θl)./C) .<= abstol)
            end
            @testset "Right tail" begin
                # set up single-grid-cell state vars
                soil
                T = [5.0]
                θl = f.(T,Tₘ,θres,θp,θw,γ) # set liquid water content according to freeze curve
                C = heatcapacity.(soil,heat,θw,θl,θm,θo)
                H = let T = T.-1.0,
                    θl = f.(T,Tₘ,θres,θp,θw,γ)
                    C = heatcapacity.(soil,heat,θw,θl,θm,θo);
                   enthalpy.(T,C,L,θl) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,dHdT=similar(C),H=H,θl=θl,)
                sfcc(soil, heat, state)
                @test all(abs.(T.-(H .- L.*θl)./C) .<= abstol)
            end
            @testset "Near zero" begin
                # set up single-grid-cell state vars
                soil
                T = [-0.05]
                θl = f.(T,Tₘ,θres,θp,θw,γ) # set liquid water content according to freeze curve
                C = heatcapacity.(soil,heat,θw,θl,θm,θo)
                H = let T = T.+0.04,
                    θl = f.(T,Tₘ,θres,θp,θw,γ)
                    C = heatcapacity.(soil,heat,θw,θl,θm,θo);
                   enthalpy.(T,C,L,θl) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,dHdT=similar(C),H=H,θl=θl,)
                sfcc(soil, heat, state)
                @test all(abs.(T.-(H .- L.*θl)./C) .<= abstol)
            end
        end
    end
    # TODO: DRY violation; a lot of this code is redundant and could possibly be
    # shared between different freeze curve tests.
    @testset "Westermann freeze curve" begin
        @testset "Sanity checks" begin
            f = Westermann() |> stripparams |> stripunits
            θtot = 0.8
            δ = 0.1
            @test isapprox(f(-10.0,Tₘ,θres,θtot,θtot,δ), 0.0, atol=1e-2)
            @test f(0.0,Tₘ,θres,θtot,θtot,δ) ≈ θtot
            θl = f(-0.1,Tₘ,θres,θtot,θtot,δ)
            @test θl > 0.0 && θl < 1.0
        end
        @testset "Newton solver checks" begin
            abstol = 1e-2
            reltol = 1e-4
            δ = 0.1
            f = Westermann() |> stripparams |> stripunits
            sfcc = SFCC(f, SFCCNewtonSolver(abstol=abstol, reltol=reltol, onfail=:error))
            heat = Heat(freezecurve=sfcc) |> stripparams |> stripunits
            L = heat.L
            @testset "Left tail" begin
                T = [-5.0]
                θl = f.(T,Tₘ,θres,θp,θw,δ) # set liquid water content according to freeze curve
                C = heatcapacity.(soil,heat,θw,θl,θm,θo)
                H = let T = T.+1.0,
                    θl = f.(T,Tₘ,θres,θp,θw,δ)
                    C = heatcapacity.(soil,heat,θw,θl,θm,θo);
                   enthalpy.(T,C,L,θl) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,dHdT=similar(C),H=H,θl=θl,)
                sfcc(soil, heat, state)
                @test all(abs.(T.-(H .- L.*θl)./C) .<= abstol)
            end
            @testset "Right tail" begin
                T = [5.0]
                θl = f.(T,Tₘ,θres,θp,θw,δ) # set liquid water content according to freeze curve
                C = heatcapacity.(soil,heat,θw,θl,θm,θo)
                H = let T = T.-1,
                    θl = f.(T,Tₘ,θres,θp,θw,δ)
                    C = heatcapacity.(soil,heat,θw,θl,θm,θo);
                   enthalpy.(T,C,L,θl) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,dHdT=similar(C),H=H,θl=θl,)
                sfcc(soil, heat, state)
                @test all(abs.(T.-(H .- L.*θl)./C) .<= abstol)
            end
            @testset "Near zero" begin
                T = [-0.05]
                θl = f.(T,Tₘ,θres,θp,θw,δ) # set liquid water content according to freeze curve
                C = heatcapacity.(soil,heat,θw,θl,θm,θo)
                H = let T = T.+0.04,
                    θl = f.(T,Tₘ,θres,θp,θw,δ)
                    C = heatcapacity.(soil,heat,θw,θl,θm,θo);
                   enthalpy.(T,C,L,θl) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,dHdT=similar(C),H=H,θl=θl,)
                sfcc(soil, heat, state)
                @test all(abs.(T.-(H .- L.*θl)./C) .<= abstol)
            end
        end
    end
    @testset "Dall'Amico freeze curve" begin
        @testset "Sanity checks" begin
            f = DallAmico() |> stripparams |> stripunits
            θsat = 0.8
            α = 4.0
            n = 2.0
            L = ustrip(Heat().L)
            @test isapprox(f(-10.0,Tₘ,θres,θsat,θsat,L,α,n), 0.0, atol=1e-6)
            @test f(0.0,Tₘ,θres,θsat,θsat,L,α,n) ≈ θsat
            θl = f(-0.1,Tₘ,θres,θsat,θsat,L,α,n)
            @test θl > 0.0 && θl < 1.0
        end
        @testset "Newton solver checks" begin
            abstol = 1e-2
            reltol = 1e-4
            θsat = 0.8
            α = 4.0
            n = 2.0
            Tₘ = 0.0
            f = DallAmico() |> stripparams |> stripunits
            sfcc = SFCC(f, SFCCNewtonSolver(abstol=abstol, reltol=reltol, onfail=:error))
            heat = Heat(freezecurve=sfcc) |> stripparams |> stripunits
            L = heat.L
            @testset "Left tail" begin
                T = [-5.0]
                θl = f.(T,Tₘ,θres,θp,θw,L,α,n) # set liquid water content according to freeze curve
                C = heatcapacity.(soil,heat,θw,θl,θm,θo)
                H = let T = T.+1.0,
                    θl = f.(T,Tₘ,θres,θp,θw,L,α,n)
                    C = heatcapacity.(soil,heat,θw,θl,θm,θo);
                   enthalpy.(T,C,L,θl) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,dHdT=similar(C),H=H,θl=θl,)
                @inferred sfcc(soil, heat, state)
                @test all(abs.(T.-(H .- L.*θl)./C) .<= abstol)
            end
            @testset "Right tail" begin
                T = [5.0]
                θl = f.(T,Tₘ,θres,θp,θw,L,α,n) # set liquid water content according to freeze curve
                C = heatcapacity.(soil,heat,θw,θl,θm,θo)
                H = let T = T.-1.0,
                    θl = f.(T,Tₘ,θres,θp,θw,L,α,n)
                    C = heatcapacity.(soil,heat,θw,θl,θm,θo);
                   enthalpy.(T,C,L,θl) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,dHdT=similar(C),H=H,θl=θl,)
                @inferred sfcc(soil, heat, state)
                @test all(abs.(T.-(H .- L.*θl)./C) .<= abstol)
            end
            @testset "Near zero" begin
                T = [-0.05]
                θl = f.(T,Tₘ,θres,θp,θw,L,α,n) # set liquid water content according to freeze curve
                C = heatcapacity.(soil,heat,θw,θl,θm,θo)
                H = let T = T.+0.04,
                    θl = f.(T,Tₘ,θres,θp,θw,L,α,n)
                    C = heatcapacity.(soil,heat,θw,θl,θm,θo);
                   enthalpy.(T,C,L,θl) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,dHdT=similar(C),H=H,θl=θl,)
                @inferred sfcc(soil, heat, state)
                @test all(abs.(T.-(H .- L.*θl)./C) .<= abstol)
            end
        end
    end
    @testset "Newton solver autodiff" begin
        # set up
        θsat = 0.8
        θres = 0.05
        Tₘ = 0.0
        γ = 0.1
        f = McKenzie() |> stripparams |> stripunits
        sfcc = SFCC(f, SFCCNewtonSolver(onfail=:error))
        heat = Heat(freezecurve=sfcc) |> stripparams |> stripunits
        L = heat.L
        T = [-0.1]
        θl = f.(T,Tₘ,θres,θp,θw,γ) # set liquid water content according to freeze curve
        C = heatcapacity.(soil,heat,θw,θl,θm,θo)
        H = enthalpy.(T.+0.09,C,L,θl) # compute enthalpy at +1 degree
        # test gradients
        p = ComponentArray(γ=γ)
        ∂f∂p = ForwardDiff.gradient(p ->  sum(f.(T,Tₘ,θres,θp,θw,p.γ)), p)
        @test all(isfinite.(∂f∂p))
        function F(p)
            T_ = similar(T,eltype(p))
            T_ .= T
            C = similar(C,eltype(p))
            θl = similar(θl,eltype(p))
            state = (T=T_,C=C,dHdT=similar(C),H=H,θl=θl,)
            @set! f.γ = p.γ
            SFCC(f)(soil, heat, state)
            state.T[1]
        end
        p = ComponentArray(γ=[γ])
        ∂F∂p = ForwardDiff.gradient(F, p)
        @test all(isfinite.(∂F∂p))
    end
end;

using BenchmarkTools
function benchmarksfcc()
    θsat = 0.8
    θres = 0.05
    α = 4.0
    n = 2.0
    Tₘ = 0.0
    f = DallAmico()
    sfcc = SFCC(f, SFCCNewtonSolver(α₀=1.0, τ=0.75, onfail=:error))
    soil = Soil()
    heat = Heat(freezecurve=sfcc) |> stripparams |> stripunits
    L = heat.L
    # set up multi-grid-cell state vars
    T = [-15.0 for i in 1:10]
    θw = Soils.soilcomponent(Val{:θw}(), soil.para)
    θp = Soils.soilcomponent(Val{:θp}(), soil.para)
    θm = Soils.soilcomponent(Val{:θm}(), soil.para)
    θo = Soils.soilcomponent(Val{:θo}(), soil.para)
    θl = f.(T,Tₘ,θres,θp,θw,L,α,n) # set liquid water content according to freeze curve
    C = heatcapacity.(soil,heat,θw,θl,θm,θo)
    H = let T = T.+14.999,
            θl = f.(T,Tₘ,θres,θp,θw,L,α,n) 
            C = heatcapacity.(soil,heat,θw,θl,θm,θo);
        enthalpy.(T,C,L,θl) # compute enthalpy at "true" temperature
    end
    state = (T=T,C=C,dHdT=similar(C),H=H,θl=θl,θw=θw.+zero(T))
    # sfcc.solver(soil, heat, state, sfcc.f, sfcc.∇f)
    # @time begin
    #     state.T .= -0.05
    #     sfcc.solver(soil, heat, state, sfcc.f, sfcc.∇f)
    # end
    result = @benchmark begin
        $sfcc($soil, $heat, $state)
    end
    show(stdout, "text/plain", result)
    println("\nsolution: $(state.T[1])")
end
