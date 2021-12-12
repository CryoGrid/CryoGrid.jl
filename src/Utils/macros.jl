"""
Similar to Unitful.@u_str (i.e. u"kg") but conditional on debug mode being enabled. Otherwise, no unit is applied.
This should be used to apply units (and thus dimensional analysis checks) to physical quantities at test time but
not during normal execution to avoid unnecessary overhead.
"""
macro xu_str(unit) CryoGrid.CRYOGRID_DEBUG ? :(@u_str($unit)) : 1 end
"""
Similar to @UT_str but produces a Float64 quantity type for the given unit if and only if debug mode is enabled.
If debug mode is not enabled, plain Float64 is used instead.
"""
macro Float_str(unit) CryoGrid.CRYOGRID_DEBUG ? :(typeof(@u_str($unit)*0.0)) : :(Float64) end
"""
Similar to @UT_str but produces a Float64 quantity type for the given unit if and only if debug mode is enabled.
If debug mode is not enabled, plain Real is used instead.
"""
macro Real_str(unit) CryoGrid.CRYOGRID_DEBUG ? :(typeof(@u_str($unit)*0.0)) : :(Real) end
"""
Similar to @UT_str but produces a Float64 quantity type for the given unit if and only if debug mode is enabled.
If debug mode is not enabled, plain Number is used instead.
"""
macro Number_str(unit) CryoGrid.CRYOGRID_DEBUG ? :(typeof(@u_str($unit)*0.0)) : :(Number) end
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