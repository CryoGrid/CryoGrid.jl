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
    θwi = Soils.soilcomponent(Val{:θwi}(), soil.para)
    θp = Soils.soilcomponent(Val{:θp}(), soil.para)
    θm = Soils.soilcomponent(Val{:θm}(), soil.para)
    θo = Soils.soilcomponent(Val{:θo}(), soil.para)
    θw = f.(T,Tₘ,θres,θp,θwi,Lf,α,n) # set liquid water content according to freeze curve
    C = heatcapacity.(soil,heat,θwi,θw,θm,θo)
    H = let T = T .+ 4.99,
            θw = f.(T,Tₘ,θres,θp,θwi,Lf,α,n),
            C = heatcapacity.(soil,heat,θwi,θw,θm,θo);
        enthalpy.(T,C,L,θw) # compute enthalpy at "true" temperature
    end
    state = (T=T,C=C,dHdT=similar(C),H=H,θw=θw,θwi=θwi)
    params = Utils.selectat(1, identity, Soils.sfccargs(f, soil, heat, state))
    f_params = tuplejoin((T[1],), params)
    # @btime $sfcc.f($f_params...)
    # @btime $sfcc.∇f($f_params...)
    # @btime Soils.residual($soil, $heat, $T[1], $H[1], $L, $sfcc.f, $params, $θwi, $θm, $θo)
    # @code_warntype HeatConduction.enthalpyinv(soil, heat, state, 1)
    # Soils.sfccsolve(solver, soil, heat, sfcc.f, sfcc.∇f, params, H[1], L, θwi, θm, θo, 0.0)
    @btime Soils.sfccsolve($solver, $soil, $heat, $sfcc.f, $sfcc.∇f, $params, $H[1], $L, $θwi, $θm, $θo, 0.0)
    # @btime Soils.residual($soil, $heat, $T[1], $H[1], $L, $sfcc.f, $params, $θwi, $θm, $θo)
    # res = @btime Soils.sfccsolve($solver, $soil, $heat, $sfcc.f, $sfcc.∇f, $params, $H[1], $L, $θwi, $θm, $θo)
    @btime freezethaw!($soil, $heat, $state)
    println("\nsolution: $(state.T[1]), $(state.θw[1])")
end
