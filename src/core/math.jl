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

export ∇,∇²,heaviside
