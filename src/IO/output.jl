"""
    CryoGridOutput{TSol}

Helper type that stores the raw output from a CryoGrid run along with `DimArray` views of all
logged variables. `CryoGridOutput` overrides `Base.getproperty` to allow for direct dot-syntax
access of state variables. For example, if your model has a grid variable named `T`, `out.T` returns a `DimArray`
with indexed time and depth axes. For OrdinaryDiffEq.jl outputs, the `ODESolution` can be accessed via `out.sol`,
or for convenience, the continuous solution at time `t` can be computed via `out(t)` which is equivalent to
`withaxes(out.sol(t))`.
"""
struct CryoGridOutput{TSol}
    ts::Vector{DateTime}
    sol::TSol
    data::NamedTuple
    CryoGridOutput(ts::Vector{DateTime}, sol::TSol, data::NamedTuple) where TSol = new{TSol}(ts, sol, data)
end
"""
Evaluates the continuous solution at time `t`.
"""
(out::CryoGridOutput)(t::Real) = CryoGrid.withaxes(out.sol(t), CryoGrid.Tile(out.sol.prob.f))
(out::CryoGridOutput)(t::DateTime) = out(Dates.datetime2epochms(t)/1000.0)
# Overrides from Base
function Base.show(io::IO, out::CryoGridOutput)
    withindent(str) = "    $str"
    countvars(x) = 1
    countvars(nt::NamedTuple) = sum(map(countvars, nt))
    format(res::Tuple) = join(res, "\n")
    describe(key, val) = withindent("$key => $(typeof(val).name.wrapper) of $(eltype(val)) with dimensions $(size(val))")
    describe(key, val::NamedTuple) = withindent("$key => \n$(format(map(withindent ∘ describe, keys(val), values(val))))")
    data = out.data
    nvars = countvars(data)
    println(io, "CryoGridOutput with $(length(out.ts)) time steps ($(out.ts[1]) to $(out.ts[end])) and $(nvars != 1 ? "$nvars variables" : "1 variable"):")
    strs = map(describe, keys(data), values(data))
    for r in strs
        println(io, r)
    end
end
DimensionalData.DimStack(out::CryoGridOutput, var::Symbol, vars::Symbol...) = DimStack((;map(n -> n => getproperty(out, n), tuple(var, vars...))...))
Base.keys(out::CryoGridOutput) = Base.propertynames(out.data)
Base.propertynames(out::CryoGridOutput) = tuple(fieldnames(typeof(out))..., propertynames(out.data)...)
function Base.getproperty(out::CryoGridOutput, sym::Symbol)
    if sym in (:sol,:ts,:data)
        getfield(out, sym)
    else
        out.data[sym]
    end
end
Base.Dict(out::CryoGridOutput) = Dict(map(k -> string(k) => getproperty(out, k), keys(out))...)

dimstr(::Ti) = "time"
dimstr(::Z) = "depth"
dimstr(::Dim{name}) where name = string(name)

function write_netcdf!(filename::String, out::CryoGridOutput, filemode="c")
    NCD.Dataset(filename, filemode) do ds
        # this assumes that the primary state variable has time and depth axes
        NCD.defDim(ds, "time", size(out.data[1], Ti))
        NCD.defDim(ds, "depth", size(out.data[1], Z))
        t = NCD.defVar(ds, "time", Float64, ("time",), attrib=Dict("units" => NCD.CFTime.DEFAULT_TIME_UNITS))
        z = NCD.defVar(ds, "depth", Float64, ("depth",))
        t[:] = collect(dims(out.data[1], Ti))
        z[:] = ustrip.(dims(out.data[1], Z))
        for var in keys(out)
            _write_ncd_var!(ds, var, getproperty(out, var))
        end
    end
end

function _write_ncd_var!(ds::NCD.Dataset, key::Symbol, data::AbstractDimArray)
    # drop dims of length 1; would be nice if there were a squeeze(..) function...
    single_dims = filter(d -> length(d) == 1, dims(data))
    data = isempty(single_dims) ? data : dropdims(data, dims=single_dims)
    datavar = NCD.defVar(ds, string(key), Float64, map(dimstr, dims(data)))
    idx = [Colon() for ax in axes(data)]
    setindex!(datavar, Array(ustrip.(data)), idx...)
end

function _write_ncd_var!(ds::NCD.Dataset, key::Symbol, nt::NamedTuple)
    for var in keys(nt)
        _write_ncd_var!(ds, Symbol("$key.$var"), nt[var])
    end
end

# Integrator saving cache

struct SaveCache{Tt}
    t::Vector{Tt}
    vals::Vector{Any}
end

function save!(cache::SaveCache, state, t)
    if length(cache.t) == 0 || cache.t[end] < t
        push!(cache.t, t)
        push!(cache.vals, state)
    end
end

function reset!(cache::SaveCache)
    resize!(cache.t, 0)
    resize!(cache.vals, 0)
end
