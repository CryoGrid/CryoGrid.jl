function DiffEqBase.step!(integrator::CGLiteIntegrator)
    cache = integrator.cache
    copyto!(cache.uprev, integrator.u)
    u = integrator.u
    du = cache.du
    t₀ = integrator.t
    p = integrator.p
    dt = integrator.dt
    t = t₀ + dt
    tile = Tiles.resolve(Tile(integrator.sol.prob.f), u, p, t)
    # explicit update, if necessary
    explicit_step!(integrator, tile, du, u, p, t)
    # implicit update for energy state
    implicit_step!(integrator, tile, du, u, p, t)
    integrator.t = t
    integrator.step += 1
    # invoke auxiliary state saving function in CryoGridProblem
    push!(tile.data.outputs.saveval, integrator.sol.prob.savefunc(tile, integrator.u, du))
    push!(tile.data.outputs.t, integrator.t)
    # save state in solution
    push!(integrator.sol.t, integrator.t)
    push!(integrator.sol.u, integrator.u)
    return nothing
end

function explicit_step!(integrator::CGLiteIntegrator, tile::Tile, du, u, p, t)
    dt = integrator.dt
    tile(du, u, p, t, dt)
    @. u += du*dt
    return u
end

function implicit_step!(integrator::CGLiteIntegrator, tile::Tile, du, u, p, t)
    # initialize local variables
    cache = integrator.cache
    ϵ = cache.resid
    H₀ = cache.uprev.H
    T_new = cache.T_new
    T_new .= zero(eltype(u))
    H = u.H
    dH = du.H
    dt = integrator.dt
    # iterative scheme for enthalpy update
    iter_count = 1
    ϵ_max = Inf
    while (ϵ_max > integrator.alg.tolerance && iter_count <= integrator.alg.maxiters) || iter_count < integrator.alg.miniters
        # compute tile state
        tile(du, u, p, t, dt)
        # extract relevant state variables
        dHdT = getvar(Val{:∂H∂T}(), tile, u; interp=false)
        Hinv = getvar(Val{:T}(), tile, u; interp=false)
        an = @view getvar(Val{:DT_an}(), tile, u; interp=false)[2:end]
        as = @view getvar(Val{:DT_as}(), tile, u; interp=false)[1:end-1]
        ap = getvar(Val{:DT_ap}(), tile, u; interp=false)
        bp = getvar(Val{:DT_bp}(), tile, u; interp=false)

        # solve linear system
        cglite_linsolve!(cache, H, Hinv, dHdT, an, as, ap, bp, dt)

        # update current state of H and dH
        @inbounds for i in eachindex(H)
            ϵ[i] = T_new[i] - Hinv[i]
            H[i] += dHdT[i]*ϵ[i]
            dH[i] = (H[i] - H₀[i]) / dt
        end

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
        integrator.alg.verbose && @warn "iteration did not converge (t = $(convert_t(t)), ϵ_max = $(maximum(abs.(ϵ))) @ $(argmax(abs.(ϵ))))"
        integrator.sol.retcode = ReturnCode.MaxIters
    end
    return nothing
end

function cglite_linsolve!(cache::LiteImplicitEulerCache, H, Hinv, dHdT, an, as, ap, bp, dt)
    # initialize variables
    H₀ = cache.uprev.H
    T_new = cache.T_new

    # compute diagonal factor and source terms
    Sp = @. cache.Sp = -dHdT / dt;
    Sc = @. cache.Sc = (H₀ - H) / dt - Sp*Hinv;

    # lienar solve --------------------------------------------------
    @. cache.A = -an;
    @. cache.B = (ap - Sp);
    @. cache.C = -as;
    @. cache.D = Sc + bp;
    Numerics.tdma_solve!(T_new, cache.A, cache.B, cache.C, cache.D);
    #---------------------------------------------------------------

    return nothing
end
