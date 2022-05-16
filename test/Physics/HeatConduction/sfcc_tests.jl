using CryoGrid
using ComponentArrays
using ForwardDiff
using Setfield
using Test

@testset "SFCC" begin
    Tₘ = 0.0
    θres = 0.0
    soil = @pstrip Soil(para=Soils.CharacteristicFractions())
    θwi = Soils.soilcomponent(Val{:θwi}(), soil.para)
    θp = Soils.soilcomponent(Val{:θp}(), soil.para)
    θm = Soils.soilcomponent(Val{:θm}(), soil.para)
    θo = Soils.soilcomponent(Val{:θo}(), soil.para)
    @testset "McKenzie freeze curve" begin
        @testset "Sanity checks" begin
            f = @pstrip McKenzie() keep_units=true
            let θsat = 0.8,
                γ = 1.0u"K",
                Tₘ = 0.0u"°C";
                @test isapprox(f(-10.0u"°C",θsat,θsat,θres,Tₘ,γ), 0.0, atol=1e-6)
                @test f(0.0u"°C",θsat,θsat,θres,Tₘ,γ) ≈ θsat
                θw = f(-0.1u"°C",θsat,θsat,θres,Tₘ,γ)
                @test θw > 0.0 && θw < 1.0
            end
        end
        @testset "Newton solver checks" begin
            abstol = 1e-2
            reltol = 1e-4
            γ = 0.1
            f = @pstrip McKenzie()
            sfcc = SFCC(f, SFCCNewtonSolver(abstol=abstol, reltol=reltol))
            heat = @pstrip Heat(freezecurve=sfcc)
            L = heat.L
            @testset "Left tail" begin
                T = [-5.0]
                θw = f.(T,θp,θwi,θres,Tₘ,γ) # set liquid water content according to freeze curve
                C = heatcapacity.(soil,heat,θwi,θw,θm,θo)
                H = let T = T.+1.0,
                    θw = f.(T,θp,θwi,θres,Tₘ,γ)
                    C = heatcapacity.(soil,heat,θwi,θw,θm,θo);
                   enthalpy.(T,C,L,θw) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,dHdT=similar(C),H=H,θw=θw,)
                freezethaw!(soil, heat, state)
                @test all(abs.(T.-(H .- L.*θw)./C) .<= abstol)
            end
            @testset "Right tail" begin
                # set up single-grid-cell state vars
                T = [5.0]
                θw = f.(T,θp,θwi,θres,Tₘ,γ) # set liquid water content according to freeze curve
                C = heatcapacity.(soil,heat,θwi,θw,θm,θo)
                H = let T = T.-1.0,
                    θw = f.(T,θp,θwi,θres,Tₘ,γ)
                    C = heatcapacity.(soil,heat,θwi,θw,θm,θo);
                   enthalpy.(T,C,L,θw) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,dHdT=similar(C),H=H,θw=θw,)
                freezethaw!(soil, heat, state)
                @test all(abs.(T.-(H .- L.*θw)./C) .<= abstol)
            end
            @testset "Near zero" begin
                # set up single-grid-cell state vars
                T = [-0.05]
                θw = f.(T,θp,θwi,θres,Tₘ,γ) # set liquid water content according to freeze curve
                C = heatcapacity.(soil,heat,θwi,θw,θm,θo)
                H = let T = T.+0.04,
                    θw = f.(T,θp,θwi,θres,Tₘ,γ)
                    C = heatcapacity.(soil,heat,θwi,θw,θm,θo);
                   enthalpy.(T,C,L,θw) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,dHdT=similar(C),H=H,θw=θw,)
                freezethaw!(soil, heat, state)
                @test all(abs.(T.-(H .- L.*θw)./C) .<= abstol)
            end
        end
    end
    # TODO: DRY violation; a lot of this code is redundant and could possibly be
    # shared between different freeze curve tests.
    @testset "Westermann freeze curve" begin
        @testset "Sanity checks" begin
            f = @pstrip Westermann() keep_units=true
            let θtot = 0.8,
                δ = 0.1u"K",
                Tₘ = 0.0u"°C";
                @test isapprox(f(-10.0u"°C",θtot,θtot,θres,Tₘ,δ), 0.0, atol=1e-2)
                @test f(0.0u"°C",θtot,θtot,θres,Tₘ,δ) ≈ θtot
                θw = f(-0.1u"°C",θtot,θtot,θres,Tₘ,δ)
                @test θw > 0.0 && θw < 1.0
            end
        end
        @testset "Newton solver checks" begin
            abstol = 1e-2
            reltol = 1e-4
            δ = 0.1
            f = @pstrip Westermann()
            sfcc = SFCC(f, SFCCNewtonSolver(abstol=abstol, reltol=reltol))
            heat = @pstrip Heat(freezecurve=sfcc)
            L = heat.L
            @testset "Left tail" begin
                T = [-5.0]
                θw = f.(T,θp,θwi,θres,Tₘ,δ) # set liquid water content according to freeze curve
                C = heatcapacity.(soil,heat,θwi,θw,θm,θo)
                H = let T = T.+1.0,
                    θw = f.(T,θp,θwi,θres,Tₘ,δ)
                    C = heatcapacity.(soil,heat,θwi,θw,θm,θo);
                   enthalpy.(T,C,L,θw) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,dHdT=similar(C),H=H,θw=θw,)
                freezethaw!(soil, heat, state)
                @test all(abs.(T.-(H .- L.*θw)./C) .<= abstol)
            end
            @testset "Right tail" begin
                T = [5.0]
                θw = f.(T,θp,θwi,θres,Tₘ,δ) # set liquid water content according to freeze curve
                C = heatcapacity.(soil,heat,θwi,θw,θm,θo)
                H = let T = T.-1,
                    θw = f.(T,θp,θwi,θres,Tₘ,δ)
                    C = heatcapacity.(soil,heat,θwi,θw,θm,θo);
                   enthalpy.(T,C,L,θw) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,dHdT=similar(C),H=H,θw=θw,)
                freezethaw!(soil, heat, state)
                @test all(abs.(T.-(H .- L.*θw)./C) .<= abstol)
            end
            @testset "Near zero" begin
                T = [-0.05]
                θw = f.(T,θp,θwi,θres,Tₘ,δ) # set liquid water content according to freeze curve
                C = heatcapacity.(soil,heat,θwi,θw,θm,θo)
                H = let T = T.+0.04,
                    θw = f.(T,θp,θwi,θres,Tₘ,δ)
                    C = heatcapacity.(soil,heat,θwi,θw,θm,θo);
                   enthalpy.(T,C,L,θw) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,dHdT=similar(C),H=H,θw=θw,)
                freezethaw!(soil, heat, state)
                @test all(abs.(T.-(H .- L.*θw)./C) .<= abstol)
            end
        end
    end
    @testset "Dall'Amico freeze curve" begin
        @testset "Sanity checks" begin
            f = @pstrip DallAmico() keep_units=true
            let θsat = 0.8,
                α = 4.0u"1/m",
                n = 2.0,
                Tₘ = 0.0u"°C";
                Lf = stripparams(Heat()).prop.Lf
                @test isapprox(f(-10.0u"°C",θsat,θsat,Lf,θres,Tₘ,α,n), 0.0, atol=1e-3)
                @test f(0.0u"°C",θsat,θsat,Lf,θres,Tₘ,α,n) ≈ θsat
                θw = f(-0.1u"°C",θsat,θsat,Lf,θres,Tₘ,α,n)
                @test θw > 0.0 && θw < 1.0
            end
        end
        @testset "Newton solver checks" begin
            abstol = 1e-2
            reltol = 1e-4
            θsat = 0.8
            α = 4.0
            n = 2.0
            Tₘ = 0.0
            f = @pstrip DallAmico()
            sfcc = SFCC(f, SFCCNewtonSolver(abstol=abstol, reltol=reltol))
            heat = @pstrip Heat(freezecurve=sfcc)
            L = heat.L
            Lf = heat.prop.Lf
            @testset "Left tail" begin
                T = [-5.0]
                θw = f.(T,θp,θwi,Lf,θres,Tₘ,α,n) # set liquid water content according to freeze curve
                C = heatcapacity.(soil,heat,θwi,θw,θm,θo)
                H = let T = T.+1.0,
                    θw = f.(T,θp,θwi,Lf,θres,Tₘ,α,n)
                    C = heatcapacity.(soil,heat,θwi,θw,θm,θo);
                   enthalpy.(T,C,L,θw) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,dHdT=similar(C),H=H,θw=θw,)
                @inferred freezethaw!(soil, heat, state)
                @test all(abs.(T.-(H .- L.*θw)./C) .<= abstol)
            end
            @testset "Right tail" begin
                T = [5.0]
                θw = f.(T,θp,θwi,Lf,θres,Tₘ,α,n) # set liquid water content according to freeze curve
                C = heatcapacity.(soil,heat,θwi,θw,θm,θo)
                H = let T = T.-1.0,
                    θw = f.(T,θp,θwi,Lf,θres,Tₘ,α,n)
                    C = heatcapacity.(soil,heat,θwi,θw,θm,θo);
                   enthalpy.(T,C,L,θw) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,dHdT=similar(C),H=H,θw=θw,)
                @inferred freezethaw!(soil, heat, state)
                @test all(abs.(T.-(H .- L.*θw)./C) .<= abstol)
            end
            @testset "Near zero" begin
                T = [-0.05]
                θw = f.(T,θp,θwi,Lf,θres,Tₘ,α,n) # set liquid water content according to freeze curve
                C = heatcapacity.(soil,heat,θwi,θw,θm,θo)
                H = let T = T.+0.04,
                    θw = f.(T,θp,θwi,Lf,θres,Tₘ,α,n)
                    C = heatcapacity.(soil,heat,θwi,θw,θm,θo);
                   enthalpy.(T,C,L,θw) # compute enthalpy at "true" temperature
                end
                state = (T=T,C=C,dHdT=similar(C),H=H,θw=θw,)
                @inferred freezethaw!(soil, heat, state)
                @test all(abs.(T.-(H .- L.*θw)./C) .<= abstol)
            end
        end
    end
    @testset "Newton solver autodiff" begin
        # set up
        θsat = 0.8
        θres = 0.05
        Tₘ = 0.0
        γ = 0.1
        f = @pstrip McKenzie()
        sfcc = SFCC(f, SFCCNewtonSolver())
        heat = @pstrip Heat(freezecurve=sfcc)
        L = heat.L
        T = [-0.1]
        θw = f.(T,θp,θwi,θres,Tₘ,γ) # set liquid water content according to freeze curve
        C = heatcapacity.(soil,heat,θwi,θw,θm,θo)
        H = enthalpy.(T.+0.09,C,L,θw) # compute enthalpy at +1 degree
        # test gradients
        p = ComponentArray(γ=γ)
        ∂f∂p = ForwardDiff.gradient(p ->  sum(f.(T,θp,θwi,θres,Tₘ,p.γ)), p)
        @test all(isfinite.(∂f∂p))
        function F(p)
            T_ = similar(T,eltype(p))
            T_ .= T
            C = similar(C,eltype(p))
            θw = similar(θw,eltype(p))
            state = (T=T_,C=C,dHdT=similar(C),H=H,θw=θw,)
            @set! heat.freezecurve.f.γ = p.γ
            freezethaw!(soil, heat, state)
            state.T[1]
        end
        p = ComponentArray(γ=[γ])
        ∂F∂p = ForwardDiff.gradient(F, p)
        @test all(isfinite.(∂F∂p))
    end
end;
