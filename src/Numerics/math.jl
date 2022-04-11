"""
    finitediff!(∂x::AbstractVector, x::AbstractVector, Δ::AbstractVector)

First order forward finite difference operator.
"""
function finitediff!(∂x::AbstractVector, x::AbstractVector, Δ::AbstractVector)
    @inbounds let x₁ = (@view x[1:end-1]),
        x₂ = (@view x[2:end]);
        @. ∂x = (x₂ - x₁) / Δ
    end
end
"""
    lineardiffusion!(∂y::AbstractVector, x::AbstractVector, Δ::AbstractVector, k::Number)

Second order Laplacian with constant diffusion k.
"""
function lineardiffusion!(∂y::AbstractVector, x::AbstractVector, Δ::AbstractVector, k::Number)
    @inbounds let x₁ = (@view x[1:end-2]),
        x₂ = (@view x[2:end-1]),
        x₃ = (@view x[3:end]),
        Δ₁ = (@view Δ[1:end-1]),
        Δ₂ = (@view Δ[2:end]);
        @. ∂y = k*((x₃ - x₂)/Δ₂ - (x₂-x₁)/Δ₁)/Δ₁
    end
end
"""
    nonlineardiffusion!(∂y, x, Δx, k, Δk)

Second order Laplacian with non-linear diffusion operator, `k`. Accelerated using `LoopVectorization.@turbo` for `Float64` vectors.
"""
function nonlineardiffusion!(∂y::AbstractVector, x::AbstractVector, Δx::AbstractVector, k::AbstractVector, Δk::AbstractVector)
    @inbounds let x₁ = (@view x[1:end-2]),
        x₂ = (@view x[2:end-1]),
        x₃ = (@view x[3:end]),
        k₁ = (@view k[1:end-1]),
        k₂ = (@view k[2:end]),
        Δx₁ = (@view Δx[1:end-1]),
        Δx₂ = (@view Δx[2:end]);
        nonlineardiffusion!(∂y, x₁, x₂, x₃, k₁, k₂, Δx₁, Δx₂, Δk)
    end
end
"""
    nonlineardiffusion(x₁, x₂, x₃, k₁, k₂, Δx₁, Δx₂, Δk)

Scalar second order Laplacian with non-linear diffusion operator, `k`.
"""
nonlineardiffusion(x₁, x₂, x₃, k₁, k₂, Δx₁, Δx₂, Δk) = (k₂*(x₃ - x₂)/Δx₂ - k₁*(x₂ - x₁)/Δx₁)/Δk
@propagate_inbounds function nonlineardiffusion!(∂y, x₁, x₂, x₃, k₁, k₂, Δx₁, Δx₂, Δk)
    @. ∂y = nonlineardiffusion(x₁, x₂, x₃, k₁, k₂, Δx₁, Δx₂, Δk)
end
@propagate_inbounds function nonlineardiffusion!(
    ∂y::AbstractVector{Float64},
    x₁::AbstractVector{Float64},
    x₂::AbstractVector{Float64},
    x₃::AbstractVector{Float64},
    k₁::AbstractVector{Float64},
    k₂::AbstractVector{Float64},
    Δx₁::AbstractVector{Float64},
    Δx₂::AbstractVector{Float64},
    Δk::AbstractVector{Float64},
)
    @turbo @. ∂y = nonlineardiffusion(x₁, x₂, x₃, k₁, k₂, Δx₁, Δx₂, Δk)
end

"""
    harmonicmean(x₁, x₂, w₁, w₂)

Simple weighted harmonic mean of two values, x₁ and x₂.
"""
harmonicmean(x₁, x₂, w₁, w₂) = (w₁ + w₂) / (w₁*x₁^-1 + w₂*x₂^-1)
"""
    harmonicmean!(h::AbstractVector, x::AbstractVector, w::AbstractVector)

Vectorized harmonic mean of elements in `x` with weights `w`. Output is stored in `h`,
which should have size `length(x)-1`.
"""
function harmonicmean!(h::AbstractVector, x::AbstractVector, w::AbstractVector)
    @inbounds let x₁ = (@view x[1:end-1]),
        x₂ = (@view x[2:end]),
        w₁ = (@view w[1:end-1]),
        w₂ = (@view w[2:end]);
        @. h = harmonicmean(x₁,x₂,w₁,w₂)
    end
end

"""
    heaviside(x)

Differentiable implementation of heaviside step function, i.e:

``h(x) = \\begin{cases} 1 & x ≥ 0 \\\\ 0 & x < 0 \\end{cases}``
"""
heaviside(x) = IfElse.ifelse(x >= zero(x), 1.0, 0.0)
"""
    logistic(x)

Numerically stable logistic function.

``σ(x) = \\begin{cases} \\frac{1}{1+\\exp(-x)} & x ≥ 0 \\\\ \\frac{\\exp(x)}{1+\\exp(x)} & x < 0 \\end{cases}``
"""
logistic(x) = IfElse.ifelse(x >= zero(x), 1 / (1 + exp(-x)), exp(x) / (1 + exp(x)))
"""
    logit(x)

Numerically stable logit function. True domain is (0,1) but inputs are
clamped to (ϵ,1-ϵ) for numerical convenience, making the effective domain
(-∞,∞).
"""
logit(x) = let x = clamp(x, eps(), 1-eps()); log(x) - log(1-x) end
"""
    softplus(x)

Numerically stable softplus function.

``s(x) = \\log(1+\\exp(-|x|)) + \\max(x,ϵ)``
"""
softplus(x) = log1p(exp(-abs(x))) + max(x,eps())
"""
    softplusinv(x)

Numerically stable softplus inverse function. True domain is (0,∞) but inputs are
clamped to (ϵ,∞) for numerical convenience, making the effective domain
(-∞,∞).
"""
softplusinv(x) = let x = clamp(x, eps(), Inf); IfElse.ifelse(x > 34, x, log(exp(x)-1)) end
# convenience functions to make composition easier, e.g:
# sqrt ∘ plusone == x -> sqrt(x. + 1.0)
minusone(x) = x .- one.(x)
plusone(x) = x .+ one.(x)

# Symbolic differentiation
"""
    ∇(f, dvar::Symbol)

Automatically generates an analytical partial derivative of `f` w.r.t `dvar` using Symbolics.jl.
To avoid symbolic tracing issues, the function should 1) be pure (no side effects or non-mathematical behavior) and 2) avoid
indeterminate control flow such as if-else or while blocks (technically should work but sometimes doesn't...). Conditional
logic can be included using `IfElse.ifelse`. Additional argument names are extracted automatically from the method signature
of `f`. Keyword arg `choosefn` should be a function which selects from available methods of `f` (returned by `methods`); defaults to `first`.
Note that `∇` uses `RuntimeGeneratedFunction` to produce a fully specialized and compiled Julia function; it may be slow on the first call
(due to compilation), but should be just as fast as handwriting it on subsequent calls.

Example:

For ``f(x,y) = 2x + xy``, ``\\frac{\\partial f}{\\partial x} = 2 + y``. Using `∇`, we can obtain this automagically:

```jldoctest
f(x,y) = 2*x + x*y
∇f_x = ∇(f,:x)
∇f_x(2.0,3.0)

# output

5.0
```
"""
function ∇(f, dvar::Symbol; choosefn=first, context_module=Numerics)
    argnames = Utils.argnames(f, choosefn)
    @assert dvar in argnames "function must have $dvar as an argument"
    dind = findfirst(s -> s == dvar, argnames)
    # Convert to symbols
    argsyms = map(s -> Symbolics.Num(SymbolicUtils.Sym{Real}(s)), argnames)
    # Generate analytical derivative of f
    x = argsyms[dind]
    ∂x = Differential(x)
    ∇f_expr = build_function(∂x(f(argsyms...)) |> expand_derivatives,argsyms...)
    ∇f = @RuntimeGeneratedFunction(context_module, ∇f_expr)
    return ∇f
end

# Function tabulation
"""
    Tabulated(f, argknots...)

Alias for `tabulate` intended for function types.
"""
Tabulated(f, argknots...) = tabulate(f, argknots...)
"""
    tabulate(f, argknots::Pair{Symbol,<:Union{Number,AbstractArray}}...)

Tabulates the given function `f` using a linear, multi-dimensional interpolant.
Knots should be given as pairs `:arg => A` where `A` is a `StepRange` or `Vector`
of input values (knots) at which to evaluate the function. `A` may also be a
`Number`, in which case a pseudo-point interpolant will be used (i.e valid on
`[A,A+ϵ]`). No extrapolation is provided by default but can be configured via
`Interpolations.extrapolate`.
"""
function tabulate(f, argknots::Pair{Symbol,<:Union{Number,AbstractArray}}...)
    initknots(a::AbstractArray) = Interpolations.deduplicate_knots!(a)
    initknots(x::Number) = initknots([x,x])
    interp(::AbstractArray) = Gridded(Linear())
    interp(::Number) = Gridded(Constant())
    extrap(::AbstractArray) = Flat()
    extrap(::Number) = Throw()
    names = map(first, argknots)
    # get knots for each argument, duplicating if only one value is provided
    knots = map(initknots, map(last, argknots))
    f_argnames = Utils.argnames(f)
    @assert all(map(name -> name ∈ names, f_argnames)) "Missing one or more arguments $f_argnames in $f"
    arggrid = Iterators.product(knots...)
    # evaluate function construct interpolant
    f = extrapolate(interpolate(Tuple(knots), map(Base.splat(f), arggrid), map(interp ∘ last, argknots)), map(extrap ∘ last, argknots))
    return f
end
function ∇(f::AbstractInterpolation)
    gradient(args...) = Interpolations.gradient(f, args...)
    return gradient
end
