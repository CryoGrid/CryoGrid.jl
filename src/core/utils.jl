"""
Similar to Unitful.@u_str (i.e. u"kg") but conditional on debug mode being enabled. Otherwise, no unit is applied.
This should be used to apply units (and thus dimensional analysis checks) to physical quantities at test time but
not during normal execution to avoid unnecessary overhead.
"""
macro xu_str(unit) CRYOGRID_DEBUG ? :(@u_str($unit)) : 1 end
"""
Similar to @UT_str but produces a Float64 quantity type for the given unit if and only if debug mode is enabled.
If debug mode is not enabled, plain Float64 is used instead.
"""
macro Float_str(unit) CRYOGRID_DEBUG ? :(typeof(@u_str($unit)*0.0)) : :(Float64) end
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
Debug ustrip. Remove units if and only if debug mode is NOT enabled.
"""
dustrip(x::Number) = CRYOGRID_DEBUG ? x : ustrip(x)
dustrip(u::Unitful.Units, x::Number) = CRYOGRID_DEBUG ? x : ustrip(u,x)

export @xu_str, @Float_str, @UFloat_str, @UT_str, dustrip

"""
Provides implementation of `Base.iterate` for structs.
"""
function structiterate(obj::A) where {A}
    names = fieldnames(A)
    if length(names) == 0; return nothing end
    gen = (getfield(obj,name) for name in names)
    (val,state) = iterate(gen)
    (val, (gen,state))
end

function structiterate(obj, state)
    gen, genstate = state
    nextitr = iterate(gen,genstate)
    isnothing(nextitr) ? nothing : (nextitr[1],(gen,nextitr[2]))
end

export structiterate

@inline tuplejoin(x) = x
@inline tuplejoin(x, y) = (x..., y...)
@inline tuplejoin(x, y, z...) = (x..., tuplejoin(y, z...)...)

export tuplejoin

const DepthAxis{D,Q} = AxisArrays.Axis{:depth,MVector{D,Q}} where {D,Q}
const ParamAxis{N} = AxisArrays.Axis{:param,MVector{N,Symbol}}
"""
    Profile{D,N,Q}

where D is the "depth" dimension or number of rows/knots defined, N is the number of parameters per row, and
Q is the type of each depth D (e.g. a quantity of meters). Profile is ultimately just a type-alias for an AxisArray
with statically defined rows and columns.
"""
const Profile{D,N,Q,T} = AxisArray{T,2,MMatrix{D,N,T},Tuple{DepthAxis{D,Q},ParamAxis{N}}} where {D,N,Q,T}
"""
    Profile(pairs...;names)

Constructs a Profile from the given pairs Q => (x1,...,xn) where x1...xn are the values defined at Q.
Column names for the resulting AxisArray can be set via the names parameter which accepts an NTuple of symbols,
where N must match the number of parameters given (i.e. n).
"""
function Profile(pairs::Pair{Q,NTuple{N,T}}...;names::Union{Nothing,NTuple{N,Symbol}}=nothing) where {N,Q,T}
    D = length(pairs)
    depths, vals = zip(pairs...)
    params = hcat(collect.(vals)...)'
    sparams = MMatrix{D,N}(params...)
    sdepths = MVector{D}(depths...)
    # auto-generate column names if not assigned
    names = isnothing(names) ? [Symbol(:x,:($i)) for i in 1:N] : names
    Profile{D,N,Q,T}(sparams,(DepthAxis{D,Q}(sdepths),ParamAxis{N}(MVector{N}(names...))))
end

"""
    interpolateprofile(profile::Profile, state; interp=Linear())

Interpolates the given profile to the corresponding variable grids. Assumes state to be indexable via the corresponding
variable symbol and that the parameter names in state and profile match.
"""
function interpolateprofile!(profile::Profile, state; interp=Linear())
    let (depths,names) = AxisArrays.axes(profile),
        z = ustrip.(depths);
        for p in names
            # in case state is unit-free, reinterpret to match eltype of profile
            pstate = reinterpret(eltype(profile),state[p])
            pgrid = state.grids[p]
            f = @> interpolate((z,), profile[:,p], Gridded(interp)) extrapolate(Flat())
            pstate .= f.(state.grids[p])   # assume length(grid) == length(state.p)
        end
    end
end

export Profile, interpolateprofile!

# Temporary fix for bug in ExprTools (see issue #14)
# TODO: Remove when fixed in official package.
function ExprTools.argument_names(m::Method)
    slot_syms = ExprTools.slot_names(m)
    arg_names = slot_syms[2:m.nargs]  # nargs includes 1 for self ref
    return arg_names
end

"""
    generate_derivative(f, dvar::Symbol)

Automatically generates an analytical partial derivative of `f` w.r.t `dvar` using ModelingToolkit/Symbolics.jl.
To avoid symbolic tracing issues, the function should 1) be pure (no side effects or non-mathematical behavior) and 2) avoid
indeterminate control flow such as if-else or while blocks (technically should work but sometimes doesn't...). Additional
argument names are extracted automatically from the method signature of `f`. Keyword arg `choosefn` should be a function
which selects from available methods of `f` (returned by `methods`); defaults to `first`.
"""
function generate_derivative(f, dvar::Symbol; choosefn=first, contextmodule=CryoGrid)
    # Parse function parameter names using ExprTools
    fms = methods(f)
    symbol(arg::Symbol) = arg
    symbol(expr::Expr) = expr.args[1]
    argnames = map(symbol, ExprTools.signature(choosefn(fms))[:args])
    @assert dvar in argnames "function must have $dvar as an argument"
    dind = findfirst(s -> s == dvar, argnames)
    # Convert to MTK symbols
    argsyms = map(s -> Num(Sym{Real}(s)), argnames)
    # Generate analytical derivative of f
    x = argsyms[dind]
    ∂x = Differential(x)
    ∇f_expr = build_function(∂x(f(argsyms...)) |> expand_derivatives,argsyms...)
    ∇f = @RuntimeGeneratedFunction(∇f_expr)
end

export generate_derivative
