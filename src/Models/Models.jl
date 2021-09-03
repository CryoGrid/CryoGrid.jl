"""
Pre-built CryoGrid configurations for rapid prototyping.
"""
module Models

using CryoGrid
using CryoGrid.InputOutput: Resource
using Statistics

include("presetgrids.jl")

Parameters = (
    # Faroux et al. doi:10.1109/IGARSS.2007.4422971
    EcoCLimMap_ULC_126_72 = Resource("EcoCLimMap_ULC_126_72", "json", "https://nextcloud.awi.de/s/nWiJr5pBoqFtw7p/download")
)

"""
    SoilLayerConfig

Helper type for representing site-specific soil layer configuration (e.g. soil and temperature profile).
"""
struct SoilLayerConfig{TSoilProfile,TTempProfile}
    soilprofile::TSoilProfile
    tempprofile::TTempProfile
end

const SamoylovDefault = SoilLayerConfig(
    # soil profile: depth => (total water, liquid water, mineral organic, porosity)
    (
        0.0u"m" => SoilComposition(xic=0.0,por=0.80,sat=1.0,org=0.75), #(θw=0.80,θm=0.05,θo=0.15,ϕ=0.80),
        0.1u"m" => SoilComposition(xic=0.0,por=0.80,sat=1.0,org=0.25), #(θw=0.80,θm=0.15,θo=0.05,ϕ=0.80),
        0.4u"m" => SoilComposition(xic=0.30,por=0.55,sat=1.0,org=0.25), #(θw=0.80,θm=0.15,θo=0.05,ϕ=0.55),
        3.0u"m" => SoilComposition(xic=0.0,por=0.50,sat=1.0,org=0.0), #(θw=0.50,θm=0.50,θo=0.0,ϕ=0.50),
        10.0u"m" => SoilComposition(xic=0.0,por=0.30,sat=1.0,org=0.0), #(θw=0.30,θm=0.70,θo=0.0,ϕ=0.30),
    ),
    TemperatureProfile(
        0.0u"m" => -1.0u"°C",
        2.0u"m" => -1.0u"°C",
        5.0u"m" => -3.0u"°C",
        10.0u"m" => -6.0u"°C",
        25.0u"m" => -9.0u"°C",
        100.0u"m" => -9.0u"°C",
        1000.0u"m" => 10.2u"°C",
    )
)

export SamoylovDefault

"""
    SoilHeat([heatvar=:H], upperbc::BoundaryProcess, soilconfig::SoilLayerConfig; grid::Grid=DefaultGrid, freezecurve::F=FreeWater()) where {F<:FreezeCurve}

Builds a simple one-layer soil/heat-conduction model with the given grid and configuration. Uses the "free water" freeze curve by default,
but this can be changed via the `freezecurve` parameter. For example, to use the Dall'Amico freeze curve, set `freezecurve=SFCC(DallAmico())`.
"""
function SoilHeat(heatvar, upperbc::BoundaryProcess, soilconfig::SoilLayerConfig;
    grid::Grid=DefaultGrid_5cm, freezecurve::F=FreeWater(), chunk_size=nothing) where {F<:FreezeCurve}
    strat = Stratigraphy(
        -2.0u"m" => top(upperbc),
        Tuple(layer[1] => subsurface(Symbol(:soil,i), Soil(comp=layer[2]), Heat(sp=heatvar,initialT=soilconfig.tempprofile, freezecurve=freezecurve)) for (i,layer) in enumerate(soilconfig.soilprofile)),
        1000.0u"m" => bottom(GeothermalHeatFlux(0.053u"J/s/m^2"))
    )
    CryoGridSetup(strat,grid,chunk_size=chunk_size)
end
SoilHeat(upperbc::BoundaryProcess, soilconfig::SoilLayerConfig; grid::Grid=DefaultGrid_2cm, freezecurve::F=FreeWater()) where {F<:FreezeCurve} =
    SoilHeat(:H, upperbc, soilconfig; grid=grid, freezecurve=freezecurve)

"""
    spinup(setup::CryoGridSetup, tspan::NTuple{2,DateTime}, p, tol, layername; kwargs...)

Implements a simple, iterative spin-up procedure.
Runs the model specified by `setup` over `tspan` until the profile mean up to `maxdepth` over the whole time span changes only within the given tolerance `tol`.
Returns the `ODESolution` generated by the final iteration.
"""
function spinup(setup::CryoGridSetup, tspan::NTuple{2,DateTime}, p, tol, layername; maxdepth=100u"m", maxiter=1000, saveat=3*3600.0, solver=Euler(), dt=60.0, solve_args...)
    prob = CryoGridProblem(setup, tspan, p)
    @info "Running initial solve ..."
    sol = solve(prob, solver, dt=dt, saveat=saveat, solve_args...)
    out = CryoGridOutput(sol)
    H = out.vars[layername].H
    grid = setup.meta[layername].grids.H
    max_ind = argmin(abs.(grid*1u"m" .- maxdepth))
    dz = Δ(grid)[1:max_ind]
    H_sub = H[1:max_ind,:]
    # transpose and dot with cell size to integrate over grid
    tspan_mean = mean(H_sub' * dz)
    tspan_mean_prev = Inf
    Δμ = abs.(tspan_mean - tspan_mean_prev)
    itercount = 1
    while Δμ > tol && itercount <= maxiter
        @info "[Iteration $itercount] tspan mean: $tspan_mean, Δμ: $Δμ"
        tspan_mean_prev = tspan_mean
        prob = remake(prob, u0=sol.u[end])
        sol = solve(prob, solver, dt=dt, saveat=saveat, solve_args...)
        out = CryoGridOutput(sol)
        H = out.vars[layername].H
        H_sub = H[1:max_ind,:]
        tspan_mean = mean(H_sub' * dz)
        Δμ = abs.(tspan_mean - tspan_mean_prev)
        itercount += 1
    end
    if Δμ > tol
        @warn "Energy state did not converge! Stopping after $maxiter iterations; tspan mean: $tspan_mean, Δμ: $Δμ"
    else
        @info "Finished after $itercount iterations; tspan mean: $tspan_mean, Δμ: $Δμ"
    end
    return out
end

end