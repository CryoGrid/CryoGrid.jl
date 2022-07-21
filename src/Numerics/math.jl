function ∇(f::F, x::Number) where {F}
    res = ForwardDiff.derivative!(ForwardDiff.DiffResult(zero(x), zero(x)), f, x)
    return res.value, res.derivs[1]
end
function ∇(f::F, x::AbstractArray) where {F}
    res = ForwardDiff.gradient!(ForwardDiff.DiffResult(eltype(x), zero(x)), f, x)
    return res.value, res.derivs
end

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
@inline nonlineardiffusion(x₁, x₂, x₃, k₁, k₂, Δx₁, Δx₂, Δk) = (k₂*(x₃ - x₂)/Δx₂ - k₁*(x₂ - x₁)/Δx₁)/Δk
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
    if USE_TURBO
        @turbo @. ∂y = nonlineardiffusion(x₁, x₂, x₃, k₁, k₂, Δx₁, Δx₂, Δk)
    else
        @. ∂y = nonlineardiffusion(x₁, x₂, x₃, k₁, k₂, Δx₁, Δx₂, Δk)
    end
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
    tdma_solve!(x, a, b, c, d)

Tridiagonal matrix solver; borrowed from CryoGridLite. Modifies all input vectors in-place.
"""
function tdma_solve!(x, a, b, c, d, factorized=false)
    #a, b, c are the column vectors for the compressed tridiagonal matrix, d is the right vector
    n = length(b); # n is the number of rows
    #x = zeros(n,1);
    if !factorized
        # Modify the first-row coefficients
        c[1] = c[1] / b[1];    # Division by zero risk.
        d[1] = d[1] / b[1];    # Division by zero would imply a singular matrix.

        @inbounds @fastmath for i = 2:n-1
            temp = b[i] - a[i] * c[i-1];
            c[i] = c[i] / temp;
            d[i] = (d[i] - a[i] * d[i-1]) / temp;
        end

        d[n] = (d[n] - a[n] * d[n-1])/( b[n] - a[n] * c[n-1]);
    end

    # Now back substitute.
    x[n] = d[n];
    @inbounds @fastmath for i = n-1:-1:1
        x[i] = d[i] - c[i] * x[i + 1];
    end

    return x
end

"""
    heaviside(x)

Differentiable implementation of heaviside step function, i.e:

``h(x) = \\begin{cases} 1 & x ≥ 0 \\\\ 0 & x < 0 \\end{cases}``
"""
heaviside(x) = IfElse.ifelse(x >= zero(x), one(x), zero(x))
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
