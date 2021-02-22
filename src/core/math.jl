"""
First order finite difference operator.
"""
function ∇(x::Tx,δ::Tδ,∂x::Tx) where {T,Tx,Tδ}
    let x₁ = (@view x[1:end-1]),
        x₂ = (@view x[2:end]);
        @. ∂x = (x₂ - x₁) / δ
    end
end

"""
Second order finite difference operator with constant diffusion k.
"""
function ∇²(x::Tx,δ::Tδ,k::T,∂x::T∂) where {T<:Real,Tx,Tk,Tδ,T∂}
    let x₁ = (@view x[1:end-2]),
        x₂ = (@view x[2:end-1]),
        x₃ = (@view x[3:end]),
        δ₁ = (@view δ[1:end-1]),
        δ₂ = (@view δ[2:end]);
        @. ∂x = k*((x₃ - x₂)/δ₂ - (x₂-x₁)/δ₁)/δ₁
    end
end

"""
Second order finite difference operator with non-constant diffusion function, k.
"""
function ∇²(x::Tx,δ::Tδ,k::Tk,∂x::T∂) where {Tx,Tk,Tδ,T∂}
    let x₁ = (@view x[1:end-2]),
        x₂ = (@view x[2:end-1]),
        x₃ = (@view x[3:end]),
        k₁ = (@view k[1:end-1]),
        k₂ = (@view k[2:end]),
        δ₁ = (@view δ[1:end-1]),
        δ₂ = (@view δ[2:end]);
        @inbounds @. ∂x = (k₂*(x₃ - x₂)/δ₂ - k₁*(x₂-x₁)/δ₁)/δ₁
    end
end

export ∇,∇²
