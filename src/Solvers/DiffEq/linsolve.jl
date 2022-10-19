import LinearSolve

struct TDMASolver <: SciMLBase.AbstractLinearAlgorithm end

LinearSolve.init_cacheval(alg::TDMASolver, A, b, u, Pl, Pr, maxiters, abstol, reltol, verbose) = Tridiagonal(similar(A))
LinearSolve.needs_concrete_A(alg::TDMASolver) = true

function SciMLBase.solve(cache::LinearSolve.LinearCache, alg::TDMASolver; kwargs...)
    factorized = !cache.isfresh
    if !factorized
        copyto!(cache.cacheval, cache.A)
    end
    a = diag(cache.cacheval, -1)
    b = diag(cache.cacheval, 0)
    c = diag(cache.cacheval, 1)
    d = cache.b
    x = cache.u
    Numerics.tdma_solve!(x, a, b, c, d, factorized)
    SciMLBase.build_linear_solution(alg,x,nothing,cache)
end
