const NumVec = AbstractArray{T,1} where {T<:Number}

"""
First order finite difference operator.
"""
function ∇(x::NumVec,Δ::NumVec,∂x::NumVec)
    let x₁ = (@view x[1:end-1]),
        x₂ = (@view x[2:end]);
        @inbounds @. ∂x = (x₂ - x₁) / Δ
    end
end

"""
Second order finite difference operator with constant diffusion k.
"""
function ∇²(x::NumVec,Δ::NumVec,k::Tk,∂y::NumVec) where {Tk<:Real}
    let x₁ = (@view x[1:end-2]),
        x₂ = (@view x[2:end-1]),
        x₃ = (@view x[3:end]),
        Δ₁ = (@view Δ[1:end-1]),
        Δ₂ = (@view Δ[2:end]);
        @inbounds @. ∂y = k*((x₃ - x₂)/Δ₂ - (x₂-x₁)/Δ₁)/Δ₁
    end
end

"""
Second order finite difference operator with non-constant diffusion function, k.
"""
function ∇²(x::NumVec,Δx::NumVec,k::NumVec,Δk::NumVec,∂y::NumVec)
    let x₁ = (@view x[1:end-2]),
        x₂ = (@view x[2:end-1]),
        x₃ = (@view x[3:end]),
        k₁ = (@view k[1:end-1]),
        k₂ = (@view k[2:end]),
        Δx₁ = (@view Δx[1:end-1]),
        Δx₂ = (@view Δx[2:end]);
        @inbounds @. ∂y = (k₂*(x₃ - x₂)/Δx₂ - k₁*(x₂-x₁)/Δx₁)/Δk
    end
end

"""
Differentiable implementation of heaviside step function.
"""
heaviside(x) = IfElse.ifelse(x >= 0.0, 1.0, 0.0)
"""
Numerically stable logistic function.
"""
logistic(x) = IfElse.ifelse(x >= 0, 1 / (1 + exp(-x)), exp(x) / (1 + exp(x)))
"""
Numerically stable logit function. True domain is (0,1) but inputs are
clamped to (ϵ,1-ϵ) for numerical convenience, making the effective domain
(-∞,∞).
"""
logit(x) = let x = clamp(x, eps(), 1-eps()); log(x) - log(1-x) end
"""
Numerically stable softplus function.
"""
softplus(x) = log1p(exp(-abs(x))) + max(x,eps())
"""
Numerically stable softplus inverse function. True domain is (0,∞) but inputs are
clamped to (ϵ,∞) for numerical convenience, making the effective domain
(-∞,∞).
"""
softplusinv(x) = let x = clamp(x, eps(), Inf); IfElse.ifelse(x > 34, x, log(exp(x)-1)) end

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
    fms = ExprTools.methods(f)
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
    return ∇f
end

export generate_derivative
