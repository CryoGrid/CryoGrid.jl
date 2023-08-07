"""
Similar to Unitful.@u_str (i.e. u"kg") but produces the type of the quantity rather than the instance. NOT conditional
on debug mode.
"""
macro UFloat_str(unit) :(typeof(@u_str($unit)*0.0)) end
"""
Similar to Unitful.@u_str (i.e. u"kg") but produces the type of the unit rather than the instance. NOT conditional
on debug mode.
"""
macro UT_str(unit) :(typeof(@u_str($unit))) end
"""
Convenience macro for setting scalar (single-element) arrays/vectors. It turns an expression of the form:
    `a.b = ...`
into
    `a.b[1] = ...`

This is primarily intended for code clarity, i.e to clearly discern scalar and non-scalar values.
"""
macro setscalar(expr)
    refexpr = expr.args[1]
    valexpr = expr.args[2]
    quote
        $(esc(refexpr))[1] = $(esc(valexpr))
    end
end
"""
Prepends `expr` with `Threads.@threads` if and only if `Threads.nthreads() > 1`, thus avoiding the overhead of
`@threads` when running in single-threaded mode.

Credit to @ranocha (Hendrik Ranocha)
https://discourse.julialang.org/t/overhead-of-threads-threads/53964/22
"""
macro threaded(expr)
    esc(quote
        if Threads.nthreads() == 1
            $(expr)
        else
            Threads.@threads $(expr)
        end
    end)
end
"""
    @pstrip(expr, kwargs...)

Convenience macro for `Utils.pstrip`; equivalent to `pstrip(expr, kwargs...)`.
"""
macro pstrip(expr, kwargs...)
    expr_esc = esc(expr)
    kwargs_esc = map(esc, kwargs)
    return :(pstrip($expr_esc; $(kwargs_esc...)))
end
"""
    sym_str(val)

Convenience macro, `sym"val"`, for creating a `Symbol` from `val`. Equivalent to `Symbol(val)`.
Use in situations where normal Julia `:val` syntax is not possible, e.g. `sym"1"` instead of `Symbol(1)`
or `sym"my.var"` instead of `Symbol("my.var")`.
"""
macro sym_str(val)
    quote
        Symbol($val)
    end
end
"""
    properties(expr)

Defines a new container type for properties which subtypes `NamedTupleWrapper`. Usage:
```julia
@properties MyProperties(
    prop1 = 1.0,
    prop2 = 2.0,
)
```
outputs:
```julia
struct MyProperties{TV} <: NamedTupleWrapper
    values::TV
    function HydraulicProperties(;
        prop1 = 1.0,
        prop2 = 2.0,
        additional_kwargs...
    )
        props = (;prop1, prop2, additional_kwargs...)
        return new{typeof(props)}(props)
    end
end
```
"""
macro properties(expr)
    @assert expr.head == :call "expression for @properties should be a method call, e.g. @properties MyProperties(a=1,b=2)"
    name, kwargs... = expr.args
    # check that all arguments are keyword arguments and throw an error otherwise
    @assert all(map(kw -> isa(kw, Expr) && kw.head == :kw, kwargs)) "non-keyword arguments are not supported by @properties"
    # get names of keyword arguments
    kwargs_names = map(kw -> esc(kw.args[1]), kwargs)
    # generate struct expression; the liberal and ugly usage of `esc` prevents the generated expression from capturing
    # variables from the declaring scope.
    quote
        struct $(esc(name)){$(esc(:TV))} <: $(esc(:(Utils.NamedTupleWrapper)))
            values::$(esc(:TV))
            function $(esc(name))(
                values::NamedTuple
            )
                new{typeof(values)}(values)
            end
            function $(esc(name))(;
                $(map(esc, kwargs)...),
                $(esc(:additional_kwargs))...
            )
                $(esc(:props)) = (; $(kwargs_names...), $(esc(:additional_kwargs))...)
                return new{$(esc(Expr(:call, :typeof, :props)))}($(esc(:props)))
            end
            $(esc(name))($(esc(:(values::NamedTuple)))) = new{$(esc(Expr(:call, :typeof, :values)))}($(esc(:values)))
        end
    end
end
