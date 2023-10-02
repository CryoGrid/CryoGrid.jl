using ..LiteImplicit

Base.@kwdef struct NLCGLite <: AbstractCryoGridNLSolverAlgorithm
    max_iter::Int = 1000
end

mutable struct NLCGLiteCache{Tu,Tt,Tcache} <: OrdinaryDiffEq.AbstractNLSolverCache
    ustep::Tu
    tstep::Tt
    innercache::Tcache
end

function build_nlcache(nlalg::NLCGLite, f, u::ComponentVector, p, t)
    innercache = LiteImplicit.LiteImplicitEulerCache(
        similar(u),
        similar(u),
        similar(u, eltype(u), length(u.H)),
        similar(u, eltype(u), length(u.H)),
        similar(u, eltype(u), length(u.H)),
        similar(u, eltype(u), length(u.H)),
        similar(u, eltype(u), length(u.H)-1),
        similar(u, eltype(u), length(u.H)),
        similar(u, eltype(u), length(u.H)-1),
        similar(u, eltype(u), length(u.H)),
    )
    return NLCGLiteCache(
        zero(u),
        zero(t),
        innercache,
    )
end

# I actually don't know what η is; it appears to be related to the residual norm
# and is probably used by the adaptive timestepping algorithms.
initial_η(nlsolver::NLSolver{<:NLCGLite}, integrator) = nlsolver.ηold

OrdinaryDiffEq.check_div(::NLCGLite) = false

OrdinaryDiffEq.@muladd function OrdinaryDiffEq.initialize!(nlsolver::NLSolver{<:NLCGLite}, integrator::DiffEqBase.DEIntegrator)
    nlsolver.cache.tstep = integrator.t + nlsolver.c * integrator.dt
    copyto!(nlsolver.cache.innercache.uprev, nlsolver.tmp)
    # this won't work if u has units
    @. nlsolver.cache.innercache.T_new = zero(eltype(integrator.u))
    @. nlsolver.cache.innercache.resid = zero(eltype(integrator.u))
    return nothing
end

OrdinaryDiffEq.@muladd function OrdinaryDiffEq.compute_step!(nlsolver::NLSolver{<:NLCGLite, true}, integrator)
    @unpack p, dt = integrator
    @unpack z, tmp, ztmp, γ, α, cache = nlsolver
    @unpack ustep, tstep, innercache = cache
    @. ustep = tmp + γ*z
    du = innercache.du
    tile = Tile(integrator)
    
    # compute tile state
    tile(du, ustep, p, tstep, dt)
    
    # update stats
    if DiffEqBase.has_stats(integrator)
        integrator.stats.nf += 1
    end
    
    # extract relevant state variables
    dHdT = getvar(Val{:∂H∂T}(), tile, ustep; interp=false)
    Hinv = getvar(Val{:T}(), tile, ustep; interp=false)
    an = @view getvar(Val{:DT_an}(), tile, ustep; interp=false)[2:end]
    as = @view getvar(Val{:DT_as}(), tile, ustep; interp=false)[1:end-1]
    ap = getvar(Val{:DT_ap}(), tile, ustep; interp=false)
    bp = getvar(Val{:DT_bp}(), tile, ustep; interp=false)
    LiteImplicit.cglite_linsolve!(innercache, ustep.H , Hinv, dHdT, an, as, ap, bp, dt)

    # update stats
    if DiffEqBase.has_stats(integrator)
        integrator.stats.nsolve += 1
    end

    # compute residual
    ϵ_max = -Inf
    @inbounds for i in eachindex(ustep.H)
        ϵ = innercache.resid[i] = innercache.T_new[i] - Hinv[i]
        dz = dHdT[i]*ϵ
        ztmp[i] = z[i] + dz
        du.H[i] = ztmp[i] / dt
        ϵ_max = max(abs(ϵ), ϵ_max)
    end

    return ϵ_max
end
