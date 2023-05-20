using OrdinaryDiffEq: NLSolver

"""
Base type for custom nonlinear solver algorithms.
"""
abstract type AbstractCryoGridNLSolverAlgorithm <: OrdinaryDiffEq.AbstractNLSolverAlgorithm end

kappa(::AbstractCryoGridNLSolverAlgorithm) = 1e-2

fast_convergence_cutoff(::AbstractCryoGridNLSolverAlgorithm) = 1//5

maxiters(::AbstractCryoGridNLSolverAlgorithm) = 100

"""
    build_nlcache(nlalg::AbstractCryoGridNLSolverAlgorithm, f, u, p, t)

Constructs the nonlinear solver cache for `nlalg`.
"""
build_nlcache(::AbstractCryoGridNLSolverAlgorithm, f, u, p, t) = error("not implemented")

function OrdinaryDiffEq.build_nlsolver(
    alg,
    nlalg::AbstractCryoGridNLSolverAlgorithm,
    u,
    uprev,
    p,
    t,
    dt,
    f,
    rate_prototype,
    ::Type{uEltypeNoUnits},
    ::Type{uBottomEltypeNoUnits},
    ::Type{tTypeNoUnits},
    γ,
    c,
    α,
    ::Val{true}
) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    # define fields of non-linear solver
    z = zero(u)
    tmp = zero(u)
    ztmp = zero(u)

    # build cache
    nlcache = build_nlcache(nlalg, f, u, p, t)

    # build non-linear solver
    ηold = one(t)
    NLSolver{true, tTypeNoUnits}(z, tmp, ztmp, γ, c, α, nlalg, kappa(nlalg),
                fast_convergence_cutoff(nlalg), ηold, 0, maxiters(nlalg),
                OrdinaryDiffEq.Divergence,
                nlcache)
end

include("nlsolve_cglite.jl")
