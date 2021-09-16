abstract type ParamTransform end
parameters(::ParamTransform) = []

struct ParamMapping{T,name,layer}
    transform::T
    ParamMapping(transform::T, name::Symbol, layer::Symbol) where {T<:ParamTransform} = new{T,name,layer}(transform)
end

struct ReparameterizedVector{T,TM,TP,M} <: DenseArray{T,1}
    pmap::TM # input/reparameterized param vector
    pout::TP # target/original param vector
    mappings::M # mapping metadata
    ReparameterizedVector(pmap::TM,pout::TP,mappings::ParamMapping...) where {T,TM<:ComponentVector{T},TP<:ComponentVector{T}} = new{T,TM,TP,typeof(mappings)}(pmap,pout,mappings)
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

function reparameterize(setup::CryoGridSetup, transforms::Pair{Symbol,<:Pair{Symbol,<:ParamTransform}}...)
    function getparam(p)
        # currently, we assume only one variable of each name in each layer;
        # this could be relaxed in the future but will need to be appropriately handled
        @assert length(p) == 1 "Found more than one parameter with name $var in $layer; this is not currently supported."
        return p[1]
    end
    model = Model(setup)
    nestedparams = map(flat(getparam; matchtype=NamedTuple), group(model, :layer, :fieldname))
    paramarr = ComponentArray(map(flat(p -> ustrip(p.val)), nestedparams))
    mappedparams = nestedparams
    mappings = ParamMapping[]
    for (layer,(var,transform)) in transforms
        @set! mappedparams[layer][var] = parameters(transform)
        push!(mappings, ParamMapping(transform, var, layer))
    end
    mappedarr = ComponentArray(map(flat(p -> ustrip(p.val)), mappedparams))
    return ReparameterizedVector(mappedarr, paramarr, mappings...)
end

@inline updateparams!(v::AbstractVector, setup::CryoGridSetup, du, u, t) = v
@inline @generated function updateparams!(rv::ReparameterizedVector{T,TM,TP,M}, setup::CryoGridSetup, du, u, t) where {T,TM,TP,M}
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
        p_out.$layer.$name = transform(p_map.$layer.$name, state, mapping.transform)
    end
end

# Transform implementations
@with_kw struct LinearTrend <: ParamTransform
    slope::Float64 = 0.0
    intercept::Float64 = 0.0
    tstart::Float64 = 0.0
    tstop::Float64 = Inf; @assert tstop > tstart
end
parameters(trend::LinearTrend) = (slope = Param(trend.slope), intercept = Param(trend.intercept))
function transform(p, state, trend::LinearTrend)
    let t = min(state.t - trend.tstart, trend.tstop),
        β = p.slope,
        α = p.intercept;
        IfElse.ifelse(t > 0, β*t + α, α)
    end
end
