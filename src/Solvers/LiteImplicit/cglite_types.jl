Base.@kwdef struct LiteImplicitEuler <: CryoGrid.CryoGridODEAlgorithm
    miniters::Int = 2
    maxiters::Int = 1000
    tolerance::Float64 = 1e-3
    verbose::Bool = true
end
# paired with CryoGridProblem OK
DiffEqBase.check_prob_alg_pairing(::CryoGridProblem, ::LiteImplicitEuler) = nothing
# error for other problems
DiffEqBase.check_prob_alg_pairing(prob, alg::LiteImplicitEuler) = throw(DiffEqBase.ProblemSolverPairingError(prob, alg))

struct LiteImplicitEulerCache{Tu,TA} <: SciMLBase.DECache
    uprev::Tu
    du::Tu
    T_new::TA
    resid::TA
    Sc::TA
    Sp::TA
    A::TA
    B::TA
    C::TA
    D::TA
end

function DiffEqBase.__init(prob::CryoGridProblem, alg::LiteImplicitEuler, args...; dt=24*3600.0, kwargs...)
    tile = Tile(prob.f)
    grid = tile.grid
    u0 = copy(prob.u0)
    t0 = prob.tspan[1]
    nsteps = Int(ceil((prob.tspan[2] - t0) / dt)) + 1
    # initialize storage
    u_storage = [copy(u0)]
    t_storage = [prob.tspan[1]]
    # evaluate tile at initial condition
    tile = Tiles.resolve(Tile(prob.f), u0, prob.p, t0)
    tile(zero(u0), u0, prob.p, t0, dt)
    # reset SavedValues on tile.data
    initialsave = prob.savefunc(tile, u0, similar(u0))
    savevals = SavedValues(Float64, typeof(initialsave))
    push!(savevals.saveval, initialsave)
    push!(savevals.t, t0)
    tile.data.outputs = savevals
    sol = CryoGridSolution(prob, u_storage, t_storage, alg, ReturnCode.Default)
    cache = LiteImplicitEulerCache(
        similar(prob.u0), # should have ComponentArray type
        similar(prob.u0),
        similar(u0, eltype(u0), length(prob.u0.H)),
        similar(u0, eltype(u0), length(prob.u0.H)),
        similar(u0, eltype(u0), length(prob.u0.H)),
        similar(u0, eltype(u0), length(prob.u0.H)),
        similar(u0, eltype(u0), length(prob.u0.H)-1),
        similar(u0, eltype(u0), length(prob.u0.H)),
        similar(u0, eltype(u0), length(prob.u0.H)-1),
        similar(u0, eltype(u0), length(prob.u0.H)),
    )
    p = isnothing(prob.p) ? prob.p : collect(prob.p)
    opts = CryoGridIntegratorOptions(dtmax=dt)
    return CryoGridIntegrator(alg, cache, opts, sol, SortedSet{typeof(t0)}(), copy(u0), p, t0, convert(eltype(prob.tspan), dt), 1, 1)
end
