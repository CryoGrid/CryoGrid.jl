struct Periodic{P,S,T} <: BoundaryProcess{P}
    period::Float"s"
    amplitude::T
    offset::T
    Periodic{P,S}(period::Q, amplitude::T=one(T), offset::T=one(T)) where
        {P<:SubSurfaceProcess,S<:BoundaryStyle,Q<:TimeQuantity,T} =
        new{P,S,T}(uconvert(u"s",period) |> dustrip, amplitude, offset)
end

(bc::Periodic)(t) = bc.amplitude*sin(π*(1/bc.period)*t) + bc.offset

BoundaryStyle(::Type{Periodic{P,S}}) where {P,S} = S()

export Periodic
