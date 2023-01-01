"""
    CryoGridEnsembleSetup{TTile<:Tile,Tkwargs}

Stores the basic model configuration for a CryoGrid single-tile ensemble.

See also [`CryoGridEnsembleProblem`](@ref), [`fitekp!`](@ref)
"""
struct CryoGridEnsembleSetup{TTile<:Tile,Tkwargs}
    tile::TTile
    tspan::NTuple{2,DateTime}
    prob_kwargs::Tkwargs
    CryoGridEnsembleSetup(tile::Tile, tspan::NTuple{2,DateTime}; prob_kwargs...) = new{typeof(tile),typeof(prob_kwargs)}(tile, tspan, prob_kwargs)
end
function default_ensemble_prob_func(setup::CryoGridEnsembleSetup, Θ::AbstractMatrix, param_map=identity)
    p = CryoGrid.parameters(setup.tile)
    u0, _ = initialcondition!(setup.tile, setup.tspan, p)
    prob = CryoGridProblem(setup.tile, u0, setup.tspan, p; setup.prob_kwargs...)
    function prob_func(prob, i, repeat)
        ϕ = param_map(Θ[:,i])
        u0, _ = initialcondition!(setup.tile, setup.tspan, ϕ)
        new_prob = remake(prob, p=ϕ, u0=u0)
        return new_prob
    end
end
function default_ensemble_output_func(outdir, prefix)
    function output_func(sol, i)
        out = CryoGridOutput(sol)
        filepath = joinpath(outdir, "$(prefix)_i=$i.jld2")
        save(filepath, out)
        return filepath
    end
    return output_func
end
"""
    CryoGridEnsembleProblem(
        setup::CryoGridEnsembleSetup,
        Θ::AbstractMatrix;
        output_dir=".",
        prob_func=default_ensemble_prob_func(setup, Θ),
        output_func=default_ensemble_output_func(output_dir, "cryogrid_ensemble_run"),
        reduction=(u,data,i) -> (append!(u,data),false),
        ensprob_kwargs...
    )

Constructs an `EnsembleProblem` from the given ensemble `setup` for a `m x N` parameter matrix `Θ`, where `N` is the size of
the ensemble and `m` is the dimensionality of the ensmeble state space (e.g. the number of parameters). `param_map` must be a
`ParameterMapping` with a transform function that accepts an `m`-dimensional vector and produces a parameter vector (or `CryoGridParams` instance)
which matches the output of `CryoGrid.parameters`. By default, `param_map` is the identity function; however, it may be customized to permit the
construction of reduced-rank or reparameterized ensembles for which the parameter space differs from the full CryoGrid model parameter space.

Keyword arguments:

`output_func`: a function `(sol,i)::Any` which processes the `ODESolution` for ensemble member `i` and returns the result.
It is recommended to save output to disk for non-trivial time spans to avoid slowdowns from serialization time when running
the ensemble using parallel workers.

`reduction`: a function `(u,data,i)` which accumulates the result of `output_func` in `u`. Defaults to just appending `data`
to `u`.

`output_dir`: Only used to specify the output directory for the default implementation of `output_func`. If a custom `output_func`
is provided, this arugment is ignored.

All additional keyword arguments will be passed to `EnsembleProblem`.

See also [`SciMLBase.EnsembleProblem`](@ref), [`CryoGridEnsembleSetup`](@ref), [`fitekp!`](@ref)
"""
function CryoGridEnsembleProblem(
    setup::CryoGridEnsembleSetup,
    Θ::AbstractMatrix;
    output_dir=".",
    prob_func=default_ensemble_prob_func(setup, Θ),
    output_func=default_ensemble_output_func(output_dir, "cryogrid_ensemble_run"),
    reduction=(u,data,i) -> (append!(u,data),false),
    ensprob_kwargs...
)
    p = CryoGrid.parameters(setup.tile)
    u0, _ = initialcondition!(setup.tile, setup.tspan, p)
    prob = CryoGridProblem(setup.tile, u0, setup.tspan, p; setup.prob_kwargs...)
    return EnsembleProblem(prob; prob_func, output_func, reduction, ensprob_kwargs...)
end
