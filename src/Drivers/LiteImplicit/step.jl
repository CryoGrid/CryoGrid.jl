function step!(integrator::CGLiteIntegrator)
    u = integrator.u
    H₀ = u.H
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
    while ϵ_max > integrator.alg.tolerance && iter_count <= integrator.alg.maxiters
        # diagnostic update and interact!
        state = getstate(tile, H, dH, t)
        fastiterate(layers(tile.strat)) do named_layer
            CryoGrid.diagnosticstep!(named_layer.obj, getproperty(state, layername(named_layer)))
        end
        stratiterate(tile.strat, state) do layer1, layer2, state1, state2
            CryoGrid.interact!(layer1, layer2, state1, state2)
        end
        k = getvar(Val{:k}(), tile, H; interp=false)
        jH = getvar(Val{:jH}(), tile, H; interp=false)
        dHdT = getvar(Val{:∂H∂T}(), tile, H; interp=false)
        Hinv = getvar(Val{:T}(), tile, H; interp=false)
        T_ub = getscalar(state.top.T_ub)

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
        bp[end] = jH[end]

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
    copyto!(integrator.u, H)
    i = integrator.step + 1
    integrator.sol.H[i] .= H
    integrator.sol.T[i] .= T_new
    integrator.sol.t[i] = t
    integrator.step = i
    integrator.t = t
    return nothing
end
