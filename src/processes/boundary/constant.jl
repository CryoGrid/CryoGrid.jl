struct Constant{P,S,T} <: BoundaryProcess{P}
    value::T
    Constant{P,S}(value::T) where {P<:SubSurfaceProcess,S<:BoundaryStyle,T} = new{P,S,T}(value)
end

# Arguments are irrelevant for Constant, so we can just use args...
(bc::Constant)(args...) = bc.value

BoundaryStyle(::Type{<:Constant{P,S}}) where {P,S} = S()

export Constant
