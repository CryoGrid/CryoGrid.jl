"""
    CryoGridOutput{TRes}

Helper type that stores the raw output from a CryoGrid run along with `DimArray` views of all
logged variables. `CryoGridOutput` overrides `Base.getproperty` to allow for direct dot-syntax
access of state variables. For example, if your model has a grid variable named `T`, `out.T` returns a `DimArray`
with indexed time and depth axes. For OrdinaryDiffEq.jl outputs, the `ODESolution` can be accessed via `out.res`,
or for convenience, the continuous solution at time `t` can be computed via `out(t)` which is equivalent to
`withaxes(out.res(t))`.
"""
struct CryoGridOutput{TRes}
    ts::Vector{DateTime}
    res::TRes
    out::NamedTuple
    CryoGridOutput(ts::Vector{DateTime}, res::TRes, out::NamedTuple) where TRes = new{TRes}(ts, res, out)
end
# Overrides from Base
function Base.show(io::IO, out::CryoGridOutput)
    countvars(x) = 1
    countvars(nt::NamedTuple) = sum(map(countvars, nt))
    repr(key, val) = "$key => $(typeof(val).name.wrapper) of $(eltype(val)) with dimensions $(size(val))"
    repr(key, val::NamedTuple) = "$key => $(map(str, keys(val), values(val)))"
    vars = out.vars
    nvars = countvars(vars)
    println(io, "CryoGridOutput with $(length(out.ts)) time steps and $(nvars != 1 ? "$nvars variables" : "1 variable"):")
    reprs = map(repr, keys(vars), values(vars))
    for r in reprs
        println(io, "   $r")
    end
end
function Base.getproperty(out::CryoGridOutput, sym::Symbol)
    if sym in (:res,:ts,:out)
        getfield(out, sym)
    else
        out.out[sym]
    end
end
