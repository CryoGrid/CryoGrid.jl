struct Constant{P,S,T} <: BoundaryProcess{P}
    value::T
    Constant(::Type{P},::Type{S}, value::T) where {P<:SubSurfaceProcess,S<:BoundaryStyle,T} = new{P,S,T}(value)
end

BoundaryStyle(::Type{Constant{P,S}}) where {P,S} = S()
