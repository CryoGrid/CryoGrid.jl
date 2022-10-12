function step!(integrator::CGLiteIntegrator)
    H₀ = integrator.u
    t₀ = integrator.t
    p = integrator.p
    dt = integrator.dt
    t = t₀ + dt
    cache = integrator.cache
    dH = cache.dH
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
    while iter_count <= integrator.alg.maxiters
        # diagnostic update and interact!
        state = getstate(tile, H, dH, t)
        fastiterate(layers(tile.strat)) do named_layer
            CryoGrid.diagnosticstep!(named_layer.obj, getproperty(state, layername(named_layer)))
        end
        stratiterate(tile.strat, state) do layer1, layer2, state1, state2
            CryoGrid.interact!(layer1, layer2, state1, state2)
        end
        k = getvar(Val{:k}(), tile, H; interp=false)
        dHdT = getvar(Val{:∂H∂T}(), tile, H; interp=false)
        Hinv = getvar(Val{:T}(), tile, H; interp=false)
        T_ub = getscalar(state.top.T_ub)
        # convergence check
        @. ϵ = T_new - Hinv
        ϵ_max = maximum(ϵ)
        if ϵ_max < integrator.alg.tolerance && iter_count >= integrator.alg.miniters
            break
        elseif !isfinite(ϵ_max)
            error("NaN values in residual at iteration $iter_count: $ϵ")
        end

        k_inner = @view k[2:end-1]
        dxpn = @view dxp[1:end-1]
        dxps = @view dxp[2:end]
        @. an = k_inner / dx / dxpn
        @. as = k_inner / dx / dxps

        # reset terms
        bp .= 0
        ap .= 0

        #Additional heat fluxes ----------------------------------------
        ap[1:end-1] .+= as
        ap[2:end] .+= an
        # account for boundary fluxes; assumes Dirichlet upper boundary and Neumann lower boundary
        ap[1] += k[1] / (dxp[1]^2 / 2)
        bp[1] = T_ub*k[1] / (dxp[1]^2 / 2)
        bp[end] = dH[end]

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
        @. H += dHdT*(T_new - Hinv)
        iter_count += 1
    end
    if iter_count > integrator.alg.maxiters && ϵ_max > integrator.alg.tolerance
        @warn "iteration did not converge (t = $(convert_t(t)), ϵ_max = $(maximum(ϵ)) @ $(argmax(ϵ)))"
    end
    # write new state into integrator
    copyto!(integrator.u, H)
    integrator.t = t
    return nothing
end
