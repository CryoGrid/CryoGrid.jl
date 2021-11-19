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
    data::NamedTuple
    CryoGridOutput(ts::Vector{DateTime}, res::TRes, data::NamedTuple) where TRes = new{TRes}(ts, res, data)
end
# Overrides from Base
function Base.show(io::IO, out::CryoGridOutput)
    countvars(x) = 1
    countvars(nt::NamedTuple) = sum(map(countvars, nt))
    describe(key, val) = "$key => $(typeof(val).name.wrapper) of $(eltype(val)) with dimensions $(size(val))"
    describe(key, val::NamedTuple) = "$key => $(map(describe, keys(val), values(val)))"
    data = out.data
    nvars = countvars(data)
    println(io, "CryoGridOutput with $(length(out.ts)) time steps and $(nvars != 1 ? "$nvars variables" : "1 variable"):")
    strs = map(describe, keys(data), values(data))
    for r in strs
        println(io, "   $r")
    end
end
function Base.getproperty(out::CryoGridOutput, sym::Symbol)
    if sym in (:res,:ts,:data)
        getfield(out, sym)
    else
        out.data[sym]
    end
end
