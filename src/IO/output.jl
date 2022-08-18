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
(out::CryoGridOutput)(t::Real) = CryoGrid.withaxes(out.sol(t), out.sol.prob.f.f)
(out::CryoGridOutput)(t::DateTime) = out(Dates.datetime2epochms(t)/1000.0)
function (out::CryoGridOutput)(t, var::Symbol)
    u = out(t)
    return getvar(var, out.sol.prob.f.f, u)
end
# Overrides from Base
function Base.show(io::IO, out::CryoGridOutput)
    countvars(x) = 1
    countvars(nt::NamedTuple) = sum(map(countvars, nt))
    format(res::Tuple) = join(res, ", ")
    describe(key, val) = "$key => $(typeof(val).name.wrapper) of $(eltype(val)) with dimensions $(size(val))"
    describe(key, val::NamedTuple) = "$key => $(format(map(describe, keys(val), values(val))))"
    data = out.data
    nvars = countvars(data)
    println(io, "CryoGridOutput with $(length(out.ts)) time steps ($(out.ts[1]) to $(out.ts[end])) and $(nvars != 1 ? "$nvars variables" : "1 variable"):")
    strs = map(describe, keys(data), values(data))
    for r in strs
        println(io, "    $r")
    end
end
stack(out::CryoGridOutput, var::Symbol, vars::Symbol...) = DimStack((;map(n -> n => getproperty(out, n), tuple(var, vars...))...))
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
