"""
    getpriors(ps::CryoGridParams)

Populates the `prior` field on all `Param`s in the `CryoGridParams` model.
"""
function getpriors(ps::CryoGridParams)
    priors = :prior ∈ keys(ps) ? ps[:prior] : Tuple(repeat([nothing], length(ps)))
    domains = :domain ∈ keys(ps) ? ps[:domain] : Tuple(repeat([nothing], length(ps)))
    return map(ps[:val], domains, priors) do value, domain, pdist
        isnothing(pdist) ? autoprior(value, one(value), upper=supremum(domain), lower=infimum(domain)) : pdist
    end
end

function CryoGridPrior(params::CryoGridParams)
    params[:prior]  = priors = getpriors(params)
    params[:param_id] = ids = param_ids(params)
    pairs = map(Pair, ids, priors)
    return PriorDistribution((; pairs...))
end
