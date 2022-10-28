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
