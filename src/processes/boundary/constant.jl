struct Constant{P,S,T} <: BoundaryProcess{P}
    value::T
    Constant{P,S}(value::T) where {P<:SubSurfaceProcess,S<:BoundaryStyle,T} = new{P,S,T}(value)
end

(bc::Constant)(t) = bc.value

BoundaryStyle(::Type{Constant{P,S}}) where {P,S} = S()

export Constant
