abstract type ParamTransform end

struct ParamMapping{T,name,layer}
    transform::T
    ParamMapping(transform::T, name::Symbol, layer::Symbol) where {T<:ParamTransform} = new{T,name,layer}(transform)
end

struct ParameterVector{T,TV,P,M} <: DenseArray{T,1}
    vals::TV # input/reparameterized param vector
    params::P # parameters grouped by layer and name
    mappings::M # mapping metadata
    ParameterVector(vals::TV,params::P,mappings::ParamMapping...) where {T,TV<:ComponentVector{T},P<:NamedTuple} = new{T,TV,P,typeof(mappings)}(vals,params,mappings)
end
vals(rv::ParameterVector) = getfield(rv, :vals)
mappings(rv::ParameterVector) = getfield(rv, :mappings)
Base.axes(rv::ParameterVector) = axes(getfield(rv, :vals))
Base.LinearIndices(rv::ParameterVector) = LinearIndices(getfield(rv, :vals))
Base.IndexStyle(::Type{<:ParameterVector}) = Base.IndexLinear()
Base.similar(rv::ParameterVector) = ParameterVector(similar(vals(rv)), similar(pout(rv)), mappings(rv))
Base.similar(rv::ParameterVector, ::Type{T}) where T = ParameterVector(similar(vals(rv), T), similar(pout(rv), T), mappings(rv))
Base.length(rv::ParameterVector) = length(getfield(rv, :vals))
Base.size(rv::ParameterVector) = size(getfield(rv, :vals))
Base.getproperty(rv::ParameterVector, sym::Symbol) = getproperty(getfield(rv, :vals), sym)
Base.getindex(rv::ParameterVector, i) = getfield(rv, :vals)[i]
Base.setproperty!(rv::ParameterVector, val, i) = setproperty!(getfield(rv, :vals), val, sym)
Base.setindex!(rv::ParameterVector, val, i) = setindex!(getfield(rv, :vals), val, i)
Base.show(io, rv::ParameterVector) = show(io, getfield(rv, :vals))
ComponentArrays.ComponentArray(rv::ParameterVector) = getfield(rv, :vals)

_paramval(p::Param) = ustrip(p.val) # extracts value from Param type and strips units
function parameterize(setup::CryoGridSetup, transforms::Pair{Symbol,<:Pair{Symbol,<:ParamTransform}}...)
    function getparam(p)
        # currently, we assume only one variable of each name in each layer;
        # this could be relaxed in the future but will need to be appropriately handled
        @assert length(p) == 1 "Found more than one parameter with name $var in $layer; this is not currently supported."
        return p[1]
    end
    model = Model(setup)
    nestedparams = mapflat(getparam, groupparams(model, :layer, :fieldname); maptype=NamedTuple)
    mappedparams = nestedparams
    mappings = ParamMapping[]
    for (layer,(var,transform)) in transforms
        @set! mappedparams[layer][var] = mapflat(getparam, groupparams(Model(transform), :fieldname); maptype=NamedTuple)
        push!(mappings, ParamMapping(transform, var, layer))
    end
    mappedarr = ComponentArray(mapflat(_paramval, mappedparams))
    return ParameterVector(mappedarr, mappedparams, mappings...)
end
@inline updateparams!(v::AbstractVector, setup::CryoGridSetup, du, u, t) = v
@inline updateparams!(v::ParameterVector{T,TV,P,Tuple{}}, setup::CryoGridSetup, du, u, t) where {T,TV,P} = v
@inline @generated function updateparams!(rv::ParameterVector{T,TV,P,M}, setup::CryoGridSetup, du, u, t) where {T,TV,P,M}
    expr = quote
        pvals = vals(rv)
        pmodel = getfield(rv, :params)
    end
    # apply parameter transforms
    for i in 1:length(M.parameters)
        push!(expr.args, :(pmodel = updateparams(pmodel, pvals, mappings(rv)[$i], setup, du, u, t)))
    end
    # flatten parameters and strip Param types
    push!(expr.args, :(return Utils.genmap(_paramval, Flatten.flatten(pmodel, ModelParameters.SELECT, ModelParameters.IGNORE))))
    return expr
end
@inline @generated function updateparams(pmodel, pvals, mapping::ParamMapping{T,name,layer}, setup::CryoGridSetup, du, u, t) where {T,name,layer}
    quote
        state = getstate(Val($(QuoteNode(layer))), setup, du, u, t)
        p = pvals.$layer.$name
        # reconstruct transform with new parameter values
        op = Flatten.reconstruct(mapping.transform, Flatten.flatten(p), ModelParameters.SELECT, ModelParameters.IGNORE)
        # apply transform and replace parameter in named tuple
        @set! pmodel.$layer.$name = Param(transform(state, op))
        return pmodel
    end
end

# Transform implementations
@with_kw struct LinearTrend{P} <: ParamTransform
    slope::P = Param(0.0)
    intercept::P = Param(0.0)
    tstart::Float64 = 0.0
    tstop::Float64 = Inf; @assert tstop > tstart
    minval::Float64 = -Inf
    maxval::Float64 = Inf
end
function transform(state, trend::LinearTrend)
    let t = min(state.t - trend.tstart, trend.tstop),
        β = trend.slope,
        α = trend.intercept;
        min(max(IfElse.ifelse(t > 0, β*t + α, α), trend.minval), trend.maxval)
    end
end
