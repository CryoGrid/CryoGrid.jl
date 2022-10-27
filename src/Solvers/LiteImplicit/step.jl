function DiffEqBase.step!(integrator::CGLiteIntegrator)
    cache = integrator.cache
    u = integrator.u
    du = cache.du
    H₀ = u.H
    t₀ = integrator.t
    p = integrator.p
    dt = integrator.dt
    t = t₀ + dt
    dH = du.H
    H = cache.H
    H .= H₀
    ϵ = cache.resid
    T_new = cache.T_new
    T_new .= zero(eltype(H))
    tile = Strat.updateparams(Tile(integrator.sol.prob.f), H, p, t)
    grid = tile.state.grid
    dx = cache.dx
    ap = cache.ap
    an = cache.an
    as = cache.as
    bp = cache.bp
    Sp = cache.Sp
    Sc = cache.Sc
    # initialize grid spacing
    dxp = Δ(grid) # grid cell thickness; length = N
    dx .= Δ(cells(grid)) # length: N - 1

    iter_count = 1
    ϵ_max = Inf
    while ϵ_max > integrator.alg.tolerance && iter_count <= integrator.alg.maxiters
        # invoke Tile step function
        CryoGrid.Strat.step!(tile, dH,  H, p, t, dt)
        dHdT = getvar(Val{:∂H∂T}(), tile, H; interp=false)
        Hinv = getvar(Val{:T}(), tile, H; interp=false)
        an = @view getvar(Val{:DT_an}(), tile, H; interp=false)[2:end]
        as = @view getvar(Val{:DT_as}(), tile, H; interp=false)[1:end-1]
        ap = getvar(Val{:DT_ap}(), tile, H; interp=false)
        bp = getvar(Val{:DT_bp}(), tile, H; interp=false)

        # bp_lat[:,j] = sum(lat_flux,dims=2)./Vp[:,j]; #[W/m³];
        # bp = bp + bp_lat[:,j];

        @. Sp = -dHdT / dt;
        @. Sc = (H₀ - H) / dt - Sp*Hinv;

        # lienar solve --------------------------------------------------
        @. cache.A = -an;
        @. cache.B = (ap - Sp);
        @. cache.C = -as;
        @. cache.D = cache.Sc + bp;
        Numerics.tdma_solve!(T_new, cache.A, cache.B, cache.C, cache.D);
        #---------------------------------------------------------------
        #update current state of H
        @. ϵ = T_new - Hinv
        @. H += dHdT*ϵ
        @. dH = (H - H₀) / dt

        # convergence check
        ϵ_max = -Inf
        for i in eachindex(ϵ)
            ϵ_max = max(ϵ_max, abs(ϵ[i]))
            if !isfinite(ϵ[i])
                error("NaN values in residual at iteration $iter_count @ t = $(convert_t(t))")
            end
        end
        iter_count += 1
    end
    if iter_count > integrator.alg.maxiters && ϵ_max > integrator.alg.tolerance
        @warn "iteration did not converge (t = $(convert_t(t)), ϵ_max = $(maximum(abs.(ϵ))) @ $(argmax(abs.(ϵ))))"
    end
    # write new state into integrator
    copyto!(integrator.u.H, H)
    i = integrator.step + 1
    push!(tile.hist.vals.saveval, integrator.sol.prob.savefunc(tile, integrator.u, get_du(integrator)))
    push!(tile.hist.vals.t, integrator.t)
    # use pre-allocated values up to time limit, then push! afterwards
    # technically step! should probably do nothing when t > tspan[2] ..?
    if i <= length(integrator.sol.u)
        integrator.sol.u[i] .= integrator.u
        integrator.sol.t[i] = t
    else
        push!(integrator.sol.u, integrator.u)
        push!(integrator.sol.t, integrator.t)
    end
    integrator.step = i
    integrator.t = t
    return nothing
end
