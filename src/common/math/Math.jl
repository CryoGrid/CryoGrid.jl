module Math

import ExprTools

using Base: @inbounds
using IfElse
using LinearAlgebra
using LoopVectorization
using RuntimeGeneratedFunctions
using Statistics
using Symbolics
using SymbolicUtils

RuntimeGeneratedFunctions.init(Math)

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
    nonlineardiffusion!(∂y::AbstractVector, x::AbstractVector, Δx::AbstractVector, k::AbstractVector, Δk::AbstractArray)

Second order Laplacian with non-linear diffusion function, k.
"""
function nonlineardiffusion!(∂y::AbstractVector, x::AbstractVector, Δx::AbstractVector, k::AbstractVector, Δk::AbstractArray)
    @inbounds let x₁ = (@view x[1:end-2]),
        x₂ = (@view x[2:end-1]),
        x₃ = (@view x[3:end]),
        k₁ = (@view k[1:end-1]),
        k₂ = (@view k[2:end]),
        Δx₁ = (@view Δx[1:end-1]),
        Δx₂ = (@view Δx[2:end]);
        @. ∂y = (k₂*(x₃ - x₂)/Δx₂ - k₁*(x₂-x₁)/Δx₁)/Δk
    end
end

"""
    nonlineardiffusion!(
        ∂y::AbstractVector{Float64},
        x::AbstractVector{Float64}, 
        Δx::AbstractVector{Float64},
        k::AbstractVector{Float64},
        Δk::AbstractVector{Float64}
    )

Second order Laplacian with non-linear diffusion function, k. Accelerated using `LoopVectorization.@turbo` for `Float64` vectors.
"""
function nonlineardiffusion!(
    ∂y::AbstractVector{Float64},
    x::AbstractVector{Float64}, 
    Δx::AbstractVector{Float64},
    k::AbstractVector{Float64},
    Δk::AbstractVector{Float64}
)
    @inbounds let x₁ = (@view x[1:end-2]),
        x₂ = (@view x[2:end-1]),
        x₃ = (@view x[3:end]),
        k₁ = (@view k[1:end-1]),
        k₂ = (@view k[2:end]),
        Δx₁ = (@view Δx[1:end-1]),
        Δx₂ = (@view Δx[2:end]);
        @turbo @. ∂y = (k₂*(x₃ - x₂)/Δx₂ - k₁*(x₂-x₁)/Δx₁)/Δk
    end
end

"""
    heaviside(x)

Differentiable implementation of heaviside step function.
"""
heaviside(x) = IfElse.ifelse(x >= 0.0, 1.0, 0.0)
"""
    logistic(x)

Numerically stable logistic function.
"""
logistic(x) = IfElse.ifelse(x >= 0, 1 / (1 + exp(-x)), exp(x) / (1 + exp(x)))
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

"""
    ∇(f, dvar::Symbol)

Automatically generates an analytical partial derivative of `f` w.r.t `dvar` using Symbolics.jl.
To avoid symbolic tracing issues, the function should 1) be pure (no side effects or non-mathematical behavior) and 2) avoid
indeterminate control flow such as if-else or while blocks (technically should work but sometimes doesn't...). Conditional
logic can be included using `IfElse.ifelse`. Additional argument names are extracted automatically from the method signature
of `f`. Keyword arg `choosefn` should be a function which selects from available methods of `f` (returned by `methods`); defaults to `first`.
"""
function ∇(f, dvar::Symbol; choosefn=first, context_module=Math)
    # Parse function parameter names using ExprTools
    fms = ExprTools.methods(f)
    symbol(arg::Symbol) = arg
    symbol(expr::Expr) = expr.args[1]
    argnames = map(symbol, ExprTools.signature(choosefn(fms))[:args])
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

export ∇

end
