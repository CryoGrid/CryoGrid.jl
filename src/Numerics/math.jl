"""
    ∇(f::F, x::Number) where {F}

Takes a function `y = f(x)` and argument `x` and returns a tuple: `(y, ∂y∂x)`.
The derivative is calculated using forward-mode automatic differentiation.
"""
function ∇(f::F, x::Number) where {F}
    res = ForwardDiff.derivative!(ForwardDiff.DiffResult(zero(x), zero(x)), f, x)
    return res.value, res.derivs[1]
end
"""
    ∇(f::F, x::AbstractArray) where {F}

Takes a function `y = f(x)` and vector-valued argument `x` and returns a tuple: `(y, ∇ₓy)`.
The gradient is calculated using forward-mode automatic differentiation.
"""
function ∇(f::F, x::AbstractArray) where {F}
    res = ForwardDiff.gradient!(ForwardDiff.DiffResult(eltype(x), zero(x)), f, x)
    return res.value, res.derivs
end
"""
    ∇(f::F) where {F}

Wraps the function `f(x)` with a method `∂f(x)` which evaluates `f` on the dual form
of `x` (i.e. converts `x` into a `ForwardDiff.Dual`) and returns the result. Note that
the derivatives are *not* extracted, so the user needs to use `ForwardDiff.value` and
`ForwardDiff.partials` to get the values and partial derivatives of all numerical
quantities produced from `x`. This method is a more flexible alternative to `∇(f, x)`
which does not assume any particular type for the output of `f`.
"""
function ∇(f::F) where {F}
    function ∂f(x; kwargs...)
        dualx = ForwardDiff.Dual{F}(x, one(x))
        return f(dualx; kwargs...)
    end
    function ∂f(x::SVector{N}; kwargs...) where {N}
        # convert each value of `x` to a ForwardDiff.Dual using `single_seed` to produce the appropriate
        # partial derivatives for each index.
        makedual(i,x) = ForwardDiff.Dual{F}(x, ForwardDiff.single_seed(ForwardDiff.Partials{N,eltype(x)}, Val{i}()))
        dualx = map(makedual, 1:N, x)
        return f(dualx; kwargs...)
    end
    return ∂f
end
# Flux calculations
@propagate_inbounds @inline _flux_kernel(x₁, x₂, Δx, k) = -k*(x₂ - x₁)/Δx
@propagate_inbounds @inline function _flux!(j::AbstractVector, x₁::AbstractVector, x₂::AbstractVector, Δx::AbstractVector, k::AbstractVector, ::Val{use_turbo}) where {use_turbo}
    @. j += _flux_kernel(x₁, x₂, Δx, k)
end
@propagate_inbounds @inline function _flux!(j::AbstractVector{Float64}, x₁::AbstractVector{Float64}, x₂::AbstractVector{Float64}, Δx::AbstractVector{Float64}, k::AbstractVector{Float64}, ::Val{true})
    @turbo @. j += _flux_kernel(x₁, x₂, Δx, k)
end
"""
    flux!(j::AbstractVector, x::AbstractVector, Δx::AbstractVector, k::AbstractVector)

Calculates the first-order, non-linear spatial flux over a discretized variable `x` with conductivity `k`.
`x` is assumed to have shape `(N,)`, `Δx` shape `(N-1,)`, and `j` and `k` shape `(N+1,)` such that `j[2:end-1]` represents
the fluxes over the inner grid cell faces. Fluxes are added to existing values in `j`.
"""
function flux!(j::AbstractVector, x::AbstractVector, Δx::AbstractVector, k::AbstractVector)
    @inbounds let x₁ = @view(x[1:end-1]),
        x₂ = @view(x[2:end]),
        k = @view(k[2:end-1]),
        j = @view(j[2:end-1]);
        _flux!(j, x₁, x₂, Δx, k, Val{USE_TURBO}())
    end
end
# Divergence
@propagate_inbounds @inline _div_kernel(j₁, j₂, Δj) = (j₁ - j₂) / Δj
@propagate_inbounds @inline function _div!(dx::AbstractVector, j₁::AbstractVector, j₂::AbstractVector, Δj::AbstractVector, ::Val{use_turbo}) where {use_turbo}
    @. dx += _div_kernel(j₁, j₂, Δj)
end
@propagate_inbounds @inline function _div!(dx::AbstractVector{Float64}, j₁::AbstractVector{Float64}, j₂::AbstractVector{Float64}, Δj::AbstractVector{Float64}, ::Val{true})
    @turbo @. dx += _div_kernel(j₁, j₂, Δj)
end
"""
    divergence!(dx::AbstractVector, j::AbstractVector, Δj::AbstractVector)

Calculates the first-order divergence over a 1D flux vector field `j` and grid cell lengths `Δj`. Divergences are added to existing values in `dx`.
"""
function divergence!(dx::AbstractVector, j::AbstractVector, Δj::AbstractVector)
    @inbounds let j₁ = @view(j[1:end-1]),
        j₂ = @view(j[2:end]);
        _div!(dx, j₁, j₂, Δj, Val{USE_TURBO}())
    end
end
# non-linear diffusion
@propagate_inbounds @inline function _nonlineardiffusion_kernel(x₁, x₂, Δx, j₁, k, Δj)
    j₂ = _flux_kernel(x₁, x₂, Δx, k)
    div = _div_kernel(j₁, j₂, Δj)
    return j₂, div
end
@propagate_inbounds @inline function _nonlineardiffusion!(dx::AbstractVector, j::AbstractVector, x::AbstractVector, Δx::AbstractVector, k::AbstractVector, Δk::AbstractVector, i::Integer)
    let x₁ = x[i-1],
        x₂ = x[i],
        j₁ = j[i-1],
        k = k[i],
        δx = Δx[i-1],
        δj = Δk[i-1];
        j₂, divx₁ = _nonlineardiffusion_kernel(x₁, x₂, δx, j₁, k, δj)
        j[i] += j₂
        dx[i-1] += divx₁
    end
end
"""
    nonlineardiffusion!(dx::AbstractVector, j::AbstractVector, x::AbstractVector, Δx::AbstractVector, k::AbstractVector, Δk::AbstractVector)

Fast alternative to `flux!` and `divergence!` which computes fluxes and divergences (via `_flux_kernel` and `_div_kernel`) in a single pass. Note, however, that
loop vectorization with `@turbo` is not possible because of necessary loop-carried dependencies. Fluxes and divergences are added to the existing values stored
in `j` and `dx`.
"""
function nonlineardiffusion!(dx::AbstractVector, j::AbstractVector, x::AbstractVector, Δx::AbstractVector, k::AbstractVector, Δk::AbstractVector) 
    @inbounds for i in 2:length(j)-1
        _nonlineardiffusion!(dx, j, x, Δx, k, Δk, i)
    end
    # compute divergence for last grid cell
    @inbounds dx[end] += _div_kernel(j[end-1], j[end], Δk[end])
    return nothing
end
# other helper functions
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
