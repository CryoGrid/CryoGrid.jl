struct Constant{P,S,T} <: BoundaryProcess{P}
    value::T
    Constant{P,S}(value::T) where {P<:SubSurfaceProcess,S<:BoundaryStyle,T} = new{P,S,T}(value)
end

# Arguments are irrelevant for Constant, so we can just use args...
(bc::Constant)(args...) = bc.value

BoundaryStyle(::Type{<:Constant{P,S}}) where {P,S} = S()

export Constant

struct Periodic{P,S,T} <: BoundaryProcess{P}
    period::Float"s"
    amplitude::T
    offset::T
    Periodic{P,S}(period::Q, amplitude::T=one(T), offset::T=one(T)) where
        {P<:SubSurfaceProcess,S<:BoundaryStyle,Q<:TimeQuantity,T} =
        new{P,S,T}(uconvert(u"s",period) |> dustrip, amplitude, offset)
end

@inline (bc::Periodic)(t) = bc.amplitude*sin(π*(1/bc.period)*t) + bc.offset
@inline (bc::Periodic)(l1,l2,p2,s1,s2) = bc(s1.t)

BoundaryStyle(::Type{<:Periodic{P,S}}) where {P,S} = S()

export Periodic
