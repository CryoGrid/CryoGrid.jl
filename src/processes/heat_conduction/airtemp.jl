struct AirTemperature{F} <: BoundaryProcess{Heat}
    forcing::F
    AirTemperature(forcing::Forcing{Float"°C"}) = new{typeof(forcing)}(forcing)
end

@inline (bc::AirTemperature)(t) = bc.forcing(t)
@inline (bc::AirTemperature)(l1,l2,p2,s1,s2) = bc(s1.t)

BoundaryStyle(::Type{<:AirTemperature}) = Dirichlet()

export AirTemperature
