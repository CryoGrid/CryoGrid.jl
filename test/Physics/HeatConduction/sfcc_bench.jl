using BenchmarkTools
using CryoGrid
using CryoGrid.Utils
using Profile
using Test

function benchmarksfcc()
    θres = 0.05
    α = 4.0
    n = 2.0
    Tₘ = 0.0
    f = DallAmico(;θres,α,n) |> stripparams |> stripunits
    solver = SFCCNewtonSolver(α₀=1.0, τ=0.75, onfail=:error)
    sfcc = SFCC(f, solver)
    soil = Soil() |> stripparams |> stripunits
    heat = Heat(freezecurve=sfcc) |> stripparams |> stripunits
    L = heat.L
    Lf = heat.prop.Lf
    # set up multi-grid-cell state vars
    T = [-5.0 for i in 1:10]
    θw = Soils.soilcomponent(Val{:θw}(), soil.para)
    θp = Soils.soilcomponent(Val{:θp}(), soil.para)
    θm = Soils.soilcomponent(Val{:θm}(), soil.para)
    θo = Soils.soilcomponent(Val{:θo}(), soil.para)
    θl = f.(T,Tₘ,θres,θp,θw,Lf,α,n) # set liquid water content according to freeze curve
    C = heatcapacity.(soil,heat,θw,θl,θm,θo)
    H = let T = T .+ 4.99,
            θl = f.(T,Tₘ,θres,θp,θw,Lf,α,n),
            C = heatcapacity.(soil,heat,θw,θl,θm,θo);
        enthalpy.(T,C,L,θl) # compute enthalpy at "true" temperature
    end
    state = (T=T,C=C,dHdT=similar(C),H=H,θl=θl,θw=θw)
    params = Utils.selectat(1, identity, Soils.sfccparams(f, soil, heat, state))
    f_params = tuplejoin((T[1],), params)
    # @btime $sfcc.f($f_params...)
    # @btime $sfcc.∇f($f_params...)
    # @btime Soils.residual($soil, $heat, $T[1], $H[1], $L, $sfcc.f, $params, $θw, $θm, $θo)
    # @code_warntype HeatConduction.enthalpyinv(soil, heat, state, 1)
    # Soils.sfccsolve(solver, soil, heat, sfcc.f, sfcc.∇f, params, H[1], L, θw, θm, θo, 0.0)
    @btime Soils.sfccsolve($solver, $soil, $heat, $sfcc.f, $sfcc.∇f, $params, $H[1], $L, $θw, $θm, $θo, 0.0)
    # @btime Soils.residual($soil, $heat, $T[1], $H[1], $L, $sfcc.f, $params, $θw, $θm, $θo)
    # res = @btime Soils.sfccsolve($solver, $soil, $heat, $sfcc.f, $sfcc.∇f, $params, $H[1], $L, $θw, $θm, $θo)
    @btime freezethaw!($soil, $heat, $state)
    println("\nsolution: $(state.T[1]), $(state.θl[1])")
end
