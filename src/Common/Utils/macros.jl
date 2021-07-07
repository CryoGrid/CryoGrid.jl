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
Similar to @UT_str but produces a Real quantity type for the given unit if and only if debug mode is enabled.
If debug mode is not enabled, plain Float64 is used instead.
"""
macro Real_str(unit) CryoGrid.CRYOGRID_DEBUG ? :(typeof(@u_str($unit)*0.0)) : :(Real) end
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

export @xu_str, @Float_str, @Real_str, @UFloat_str, @UT_str, @setscalar