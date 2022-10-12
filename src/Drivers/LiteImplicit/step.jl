function step!(integrator::CGLiteIntegrator)
    H₀ = integrator.u
    t₀ = integrator.t
    p = integrator.p
    dt = integrator.dt
    t = t₀ + dt
    cache = integrator.cache
    ϵ = cache.resid
    T_new = cache.T_new
    dH = cache.dH
    H = cache.H
    H .= H₀
    tile = prob.tile
    grid = tile.grid
    dxns = cache.dxns
    ap = cache.ap
    ans = cache.ans
    bp = cache.bp
    # initialize grid spacing
    dxp = Δ(grid) # grid cell thickness; length = N
    @. dxns = Δ(cells(grid)) # length: N - 1

    iter_count = 1
    while integrator.miniters <= iter_count <= integrator.maxiters
        # diagnostic update and interact!
        state = getstate(tile, H, dH, t)
        # compute diagnostic step
        fastiterate(layers(strat)) do named_layer
            CryoGrid.diagnosticstep!(named_layer.obj, getproperty(state, layername(named_layer)))
        end
        k = getvar(Val{:k}(), tile, H; interp=false)
        dHdT = getvar(Val{:∂H∂T}(), tile, H; interp=false)
        Hinv = getvar(Val{:T}(), tile, H; interp=false)
        T_ub = getscalar(state.top.T_ub)
        # convergence check
        @. ϵ = T_new - Hinv
        ϵ_max = maximum(ϵ)
        if ϵ_max < integrator.tolerance
            break
        end

        @. ans = k[2:end-1] / dxns / dxp;

        # reset terms
        bp .= 0
        ap .= 0

        #Additional heat fluxes----------------------------------------
        @. ap[1:end-1] += ans
        @. ap[2:end] += ans
        # account for boundary fluxes; assumes Dirichlet upper boundary and Neumann lower boundary
        ap[1] += k[1] / (dxp[1]^2 / 2)
        bp[1] = T_ub*k[1] / (dxp[1]^2 / 2)
        bp[end] = dH[end]

        # bp_lat[:,j] = sum(lat_flux,dims=2)./Vp[:,j]; #[W/m³];
        # bp = bp + bp_lat[:,j];

        @. cache.Sp = -dHdT / dt;
        @. cache.Sc = (H₀ - H) / dt - Sp*Hinv;

        # lienar solve --------------------------------------------------
        @. cache.A = -ans;
        @. cache.B = (ap - Sp);
        @. cache.C = -ans;
        @. cache.D = cache.Sc + bp;
        Numerics.tdma_solve!(T_new, cache.A, cache.B, cache.C, cache.D);
        #---------------------------------------------------------------
        #update current state of H
        H += dHdT.*(T_new - Hinv);
        iter_count += 1
    end
    if iter_count > integrator.maxiters && ϵ_max > integrator.tolerance
        @warn "iteration did not converge (t = $(convert_t(t)), ϵ_max = $(maximum(ϵ)) @ $(argmax(ϵ)))"
    end
    # write new state into integrator
    copyto!(integrator.u, H)
    integrator.t = t
end
