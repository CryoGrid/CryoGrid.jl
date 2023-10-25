"""
    CryoGridParameterEnsemble(
        prob::CryoGridProblem,
        Θ::AbstractMatrix;
        transform=identity,
        output_func=(sol, i) -> CryoGridOutput(sol),
        reduction=(u,data,i) -> (append!(u,data),false),
        ensprob_kwargs...
    )

Constructs an `EnsembleProblem` for a `m x N` parameter matrix `Θ`, where `N` is the size of the ensemble and `m` is the dimensionality of the
ensmeble parameter space (e.g. the number of parameters). `transform` should be a transform function that accepts an `m`-dimensional vector and
produces a parameter vector (or `CryoGridParams` instance) which matches the output of `CryoGrid.parameters`. By default, `param_map` is the
identity function; however, it may be customized to permit the construction of reduced-rank or reparameterized ensembles for which the parameter
space differs from the full CryoGrid model parameter space.

Keyword arguments:

`output_func`: a function `(sol,i)::Any` which processes the `ODESolution` for ensemble member `i` and returns the result.
It is recommended to save output to disk for non-trivial time spans to avoid slowdowns from serialization time when running
the ensemble using parallel workers.

`reduction`: a function `(u,data,i)` which accumulates the result of `output_func` in `u`. Defaults to just appending `data`
to `u`.

`output_dir`: Only used to specify the output directory for the default implementation of `output_func`. If a custom `output_func`
is provided, this arugment is ignored.

All additional keyword arguments will be passed to `EnsembleProblem`.

See also [`SciMLBase.EnsembleProblem`](@ref)
"""
function CryoGridParameterEnsemble(
    prob::CryoGridProblem,
    Θ::AbstractMatrix;
    transform=identity,
    output_func=(sol, i) -> CryoGridOutput(sol),
    reduction=(u,data,i) -> (append!(u,data),false),
    ensprob_kwargs...
)
    prob_func(prob::CryoGridProblem, i, repeat) = remake(prob, p=transform(Θ[:,i]))
    p = CryoGrid.parameters(prob.f)
    @assert length(p) == size(Θ,1) "The number of rows in θ ($(size(Θ,1))) does not match the number of model parameters $(length(p))."
    return EnsembleProblem(prob; prob_func, output_func, reduction, ensprob_kwargs...)
end
