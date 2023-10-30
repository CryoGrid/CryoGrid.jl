"""
    CGEuler <: SciMLBase.AbstractODEAlgorithm

Simple, lightweight implementation of the forward Euler integration algorithm.
Does not include support for fancier features such as interpolation, adaptive timestepping,
or event handling. In order to get these features, you must install the `OrdinaryDiffEq`
package.
"""
struct CGEuler <: CryoGridODEAlgorithm end

DiffEqBase.check_prob_alg_pairing(::CryoGridProblem, ::CGEuler) = nothing
DiffEqBase.check_prob_alg_pairing(prob, alg::CGEuler) = throw(DiffEqBase.ProblemSolverPairingError(prob, alg))

# if no solver is specified, solve with CGEuler
DiffEqBase.solve(prob::CryoGridProblem; kwargs...) = solve(prob, CGEuler(); kwargs...)

struct CGEulerCache{Tu} <: SciMLBase.DECache
    uprev::Tu
    du::Tu
end

function DiffEqBase.__init(prob::CryoGridProblem, alg::CGEuler, args...; dt=60.0, kwargs...)
    tile = Tile(prob.f)
    u0 = copy(prob.u0)
    du0 = zero(u0)
    t0 = prob.tspan[1]
    # initialize storage
    u_storage = [copy(u0)]
    t_storage = [prob.tspan[1]*one(eltype(u0))]
    # evaluate tile at initial condition
    tile = Tiles.resolve(Tile(prob.f), prob.p, t0)
    tile(du0, u0, prob.p, t0, dt)
    # reset SavedValues on tile.data
    initialsave = prob.savefunc(tile, u0, similar(u0))
    savevals = SavedValues(Float64, typeof(initialsave))
    push!(savevals.saveval, initialsave)
    push!(savevals.t, t0)
    tile.data.outputs = savevals
    sol = CryoGridSolution(prob, u_storage, t_storage, alg, ReturnCode.Default)
    cache = CGEulerCache(
        similar(prob.u0),
        similar(prob.u0),
    )
    p = isnothing(prob.p) ? prob.p : collect(prob.p)
    opts = CryoGridIntegratorOptions(;kwargs...)
    return CryoGridIntegrator(alg, cache, opts, sol, copy(u0), p, t0*one(eltype(u0)), dt*one(eltype(u0)), 1, 1)
end

function perform_step!(integrator::CryoGridIntegrator{CGEuler})
    copyto!(integrator.cache.uprev, integrator.u)
    u = integrator.u
    du = integrator.cache.du
    t₀ = integrator.t
    p = integrator.p
    initial_tile = Tile(integrator.sol.prob.f)
    tile = Tiles.resolve(initial_tile, p, t₀)
    # compute time derivative du
    tile(du, u, p, t₀)
    dt = integrator.dt
    # update u
    @inbounds @. u += dt*du
    integrator.u = u
    integrator.t = t₀ + dt
    integrator.dt = dt
    integrator.step += 1
    # set next dt
    # compute maximum timestep
    dtmax = CryoGrid.timestep(tile, du, u, p, t₀)
    integrator.dt = dtmax
    return nothing
end
