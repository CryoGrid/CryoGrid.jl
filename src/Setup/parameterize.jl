abstract type ParamTransform end

struct ParamMapping{T,name,layer}
    transform::T
    ParamMapping(transform::T, name::Symbol, layer::Symbol) where {T<:ParamTransform} = new{T,name,layer}(transform)
end

struct ReparameterizedVector{T,TM,TP,P,M} <: DenseArray{T,1}
    pmap::TM # input/reparameterized param vector
    pout::TP # target/original param vector
    params::P # parameters grouped by layer and name
    mappings::M # mapping metadata
    ReparameterizedVector(pmap::TM,pout::TP,params::P,mappings::ParamMapping...) where {T,TM<:ComponentVector{T},TP<:ComponentVector{T},P<:NamedTuple} = new{T,TM,TP,P,typeof(mappings)}(pmap,pout,params,mappings)
end
pmap(rv::ReparameterizedVector) = getfield(rv, :pmap)
pout(rv::ReparameterizedVector) = getfield(rv, :pout)
mappings(rv::ReparameterizedVector) = getfield(rv, :mappings)
Base.axes(rv::ReparameterizedVector) = axes(getfield(rv, :pmap))
Base.LinearIndices(rv::ReparameterizedVector) = LinearIndices(getfield(rv, :pmap))
Base.IndexStyle(::Type{<:ReparameterizedVector}) = Base.IndexLinear()
Base.similar(rv::ReparameterizedVector) = ReparameterizedVector(similar(pmap(rv)), similar(pout(rv)), mappings(rv))
Base.similar(rv::ReparameterizedVector, ::Type{T}) where T = ReparameterizedVector(similar(pmap(rv), T), similar(pout(rv), T), mappings(rv))
Base.length(rv::ReparameterizedVector) = length(getfield(rv, :pmap))
Base.size(rv::ReparameterizedVector) = size(getfield(rv, :pmap))
Base.getproperty(rv::ReparameterizedVector, sym::Symbol) = getproperty(getfield(rv, :pmap), sym)
Base.getindex(rv::ReparameterizedVector, i) = getfield(rv, :pmap)[i]
Base.setproperty!(rv::ReparameterizedVector, val, i) = setproperty!(getfield(rv, :pmap), val, sym)
Base.setindex!(rv::ReparameterizedVector, val, i) = setindex!(getfield(rv, :pmap), val, i)
Base.show(io, rv::ReparameterizedVector) = show(io, getfield(rv, :pmap))
ComponentArrays.ComponentArray(rv::ReparameterizedVector) = getfield(rv, :pmap)

function parameterize(setup::CryoGridSetup, transforms::Pair{Symbol,<:Pair{Symbol,<:ParamTransform}}...)
    function getparam(p)
        # currently, we assume only one variable of each name in each layer;
        # this could be relaxed in the future but will need to be appropriately handled
        @assert length(p) == 1 "Found more than one parameter with name $var in $layer; this is not currently supported."
        return p[1]
    end
    model = Model(setup)
    nestedparams = map(flat(getparam; matchtype=NamedTuple), group(model, :layer, :fieldname))
    mappedparams = nestedparams
    mappings = ParamMapping[]
    for (layer,(var,transform)) in transforms
        @set! mappedparams[layer][var] = group(Model(transform), :fieldname)
        push!(mappings, ParamMapping(transform, var, layer))
    end
    outarr = ComponentArray(map(flat(p -> ustrip(p.val)), nestedparams))
    mappedarr = ComponentArray(map(flat(p -> ustrip(p.val)), mappedparams))
    return ReparameterizedVector(mappedarr, outarr, mappedparams, mappings...)
end

@inline updateparams!(v::AbstractVector, setup::CryoGridSetup, du, u, t) = v
@inline @generated function updateparams!(rv::ReparameterizedVector{T,TM,TP,P,M}, setup::CryoGridSetup, du, u, t) where {T,TM,TP,P,M}
    expr = Expr(:block)
    for i in 1:length(M.parameters)
        push!(expr.args, :(updateparams!(rv, mappings(rv)[$i], setup, du, u, t)))
    end
    push!(expr.args, :(return pout(rv)))
    return expr
end
# TODO: Roll this function into the one above. This forces the compiler to compile separate functions for each individual param mapping (i.e. every parameter)
# which incurs a pretty hefty compile time cost.
@inline @generated function updateparams!(rv::ReparameterizedVector, mapping::ParamMapping{T,name,layer}, setup::CryoGridSetup, du, u, t) where {T,name,layer}
    quote
        p_map = pmap(rv)
        p_out = pout(rv)
        state = getstate(Val($(QuoteNode(layer))), setup, du, u, t)
        op = Flatten.reconstruct(mapping.transform, p_map.$layer.$name, ModelParameters.SELECT, ModelParameters.IGNORE)
        p_out.$layer.$name = transform(state, op)
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
