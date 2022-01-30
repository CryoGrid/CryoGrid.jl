abstract type ParamTransform end
shortname(::ParamTransform) = :transform
# mapping
struct ParamMapping{T,name,layer}
    transform::T
    ParamMapping(transform::T, name::Symbol, layer::Symbol) where {T<:ParamTransform} = new{T,name,layer}(transform)
end
# parameter array
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
Base.similar(rv::ParameterVector) = ParameterVector(similar(vals(rv)), getfield(rv, :params), mappings(rv)...)
Base.similar(rv::ParameterVector, ::Type{T}) where T = ParameterVector(similar(vals(rv), T), getfield(rv, :params), mappings(rv)...)
Base.length(rv::ParameterVector) = length(getfield(rv, :vals))
Base.size(rv::ParameterVector) = size(getfield(rv, :vals))
Base.getproperty(rv::ParameterVector, sym::Symbol) = getproperty(getfield(rv, :vals), sym)
Base.getindex(rv::ParameterVector, i) = getfield(rv, :vals)[i]
Base.setproperty!(rv::ParameterVector, val, i) = setproperty!(getfield(rv, :vals), val, sym)
Base.setindex!(rv::ParameterVector, val, i) = setindex!(getfield(rv, :vals), val, i)
Base.show(io::IO, ::MIME"text/plain", rv::ParameterVector{T,<:Any,<:Any,Tuple{}}) where {T} = println(io, "$(length(rv))-element ParameterVector{T}:\n$(getfield(rv, :vals))")
Base.show(io::IO, ::MIME"text/plain", rv::ParameterVector{T,<:Any,<:Any}) where {T} = println(io, "$(length(rv))-element ParameterVector{T} with $(length(mappings(rv))) mappings\n$(mappings(rv)):\n$(getfield(rv, :vals))")
ComponentArrays.ComponentArray(rv::ParameterVector) = getfield(rv, :vals)

_paramval(x) = ustrip(x)
_paramval(x::Param) = ustrip(x.val)
_paramval(x::Tuple) = collect(_paramval.(x))
function parameters(tile::Tile, transforms::Pair{Symbol,<:Pair{Symbol,<:ParamTransform}}...)
    getparam(x) = x
    getparam(x::Union{<:AbstractVector,<:Tuple}) = length(x) == 1 ? getparam(x[1]) : Tuple(getparam.(x))
    m = Model(tile)
    nestedparams = mapflat(getparam, groupparams(m, :layer, :fieldname); maptype=NamedTuple)
    mappedparams = nestedparams
    mappings = ParamMapping[]
    for (layer,(var,transform)) in transforms
        m_transform = Model(transform)
        m_transform[:transform] = repeat([shortname(transform)], length(ModelParameters.params(m_transform)))
        @set! mappedparams[layer][var] = mapflat(getparam, groupparams(m_transform, :transform, :fieldname); maptype=NamedTuple)
        push!(mappings, ParamMapping(transform, var, layer))
    end
    mappedarr = ComponentArray(mapflat(_paramval, mappedparams; maptype=NamedTuple))
    return ParameterVector(mappedarr, mappedparams, mappings...)
end
@inline @generated function updateparams!(v::AbstractVector, tile::Tile, u, du, t)
    quote
        p = ModelParameters.update(ModelParameters.params(tile), v)
        return Utils.genmap(_paramval, p)
    end
end
@inline @generated function updateparams!(rv::ParameterVector{T,TV,P,M}, tile::Tile, u, du, t) where {T,TV,P,M}
    expr = quote
        pvals = vals(rv)
        pmodel = ModelParameters.update(getfield(rv, :params), pvals)
    end
    # apply parameter transforms
    for i in 1:length(M.parameters)
        push!(expr.args, :(pmodel = _updateparam(pmodel, mappings(rv)[$i], tile, u, du, t)))
    end
    # flatten parameters and strip Param types
    push!(expr.args, :(return Utils.genmap(_paramval, ModelParameters.params(pmodel))))
    return expr
end
@inline @generated function _updateparam(pmodel, mapping::ParamMapping{T,name,layer}, tile::Tile, u, du, t) where {T,name,layer}
    quote
        state = getstate(Val($(QuoteNode(layer))), tile, u, du, t)
        p = pmodel.$layer.$name
        # reconstruct transform with new parameter values
        op = ModelParameters.update(mapping.transform, ModelParameters.params(p))
        # apply transform and replace parameter in named tuple
        @set! pmodel.$layer.$name = Param(transform(state, op))
        return pmodel
    end
end
# Transform implementations
"""
    LinearTrend{P} <: ParamTransform

Applies a linear trend to a parameter `p` by reparameterizing it as: `p = p₁*t + p₀`
"""
@with_kw struct LinearTrend{P} <: ParamTransform
    slope::P = Param(0.0)
    intercept::P = Param(0.0)
    tstart::Float64 = 0.0
    tstop::Float64 = Inf; @assert tstop > tstart
    period::Float64 = 1.0; @assert period > 0.0
    minval::Float64 = -Inf
    maxval::Float64 = Inf
end
shortname(::LinearTrend) = :trend
function transform(state, trend::LinearTrend)
    let t = min(state.t - trend.tstart, trend.tstop),
        β = trend.slope / trend.period,
        α = trend.intercept;
        min(max(β*t + α, trend.minval), trend.maxval)
    end
end
"""
Helper type for PiecewiseLinear.
"""
struct PiecewiseKnot{Tw,Tv}
    binwidth::Tw
    value::Tv
end
"""
    PiecewiseLinear{N,Tw,Tv} <: ParamTransform

Reparameterizes parameter `p` as `p = p₁δ₁t + ⋯ + pₖδₖt` where δₖ are indicators
for when `tₖ₋₁ <= t <= tₖ`. To facilitate sampling and optimization, change points
tᵢ are parameterized as bin widths, which should be strictly positive. `PiecewiseLinear`
will normalize them and scale by the size of the time interval.
"""
@with_kw struct PiecewiseLinear{N,Tw,Tv} <: ParamTransform
    initialvalue::Tv
    knots::NTuple{N,PiecewiseKnot{Tw,Tv}}
    tstart::Float64 = 0.0; @assert isfinite(tstart)
    tstop::Float64 = 1.0; @assert tstop > tstart; @assert isfinite(tstop)
end
PiecewiseLinear(initialvalue::T; tstart=0.0, tstop=1.0) where T = PiecewiseLinear{0,Nothing,T}(initialvalue, (), tstart, tstop)
PiecewiseLinear(initialvalue, knots...; tstart=0.0, tstop=1.0) = PiecewiseLinear(initialvalue, Tuple(map(Base.splat(PiecewiseKnot), knots)), tstart, tstop)
shortname(::PiecewiseLinear) = :linear
function transform(state, pc::PiecewiseLinear)
    function binindex(values::Tuple, st, en, x)
        mid = Int(floor((st + en)/2))
        if values[mid] >= x && en - st <= 1
            return st
        elseif values[mid] < x && en - st <= 1
            return en-1
        elseif values[mid] >= x
            return binindex(values, st, mid, x)
        else
            return binindex(values, mid, en, x)
        end
    end
    let tspan = pc.tstop - pc.tstart,
        t = min(max(state.t - pc.tstart, zero(state.t)), tspan),
        bins = map(k -> k.binwidth, pc.knots),
        vals = (pc.initialvalue, map(k -> k.value, pc.knots)...),
        ts = (0.0, cumsum((bins ./ sum(bins)).*tspan)...),
        i = binindex(ts, 1, length(ts), t);
        vals[i] + (vals[i+1] - vals[i])*(t - ts[i]) / (ts[i+1] - ts[i])
    end
end
