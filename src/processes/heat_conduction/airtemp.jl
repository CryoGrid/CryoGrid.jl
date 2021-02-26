struct AirTemperature{F} <: BoundaryProcess{Heat}
    forcing::F
    AirTemperature(forcing::Forcing{Float"°C"}) = new{typeof(forcing)}(forcing)
end

(bc::AirTemperature)(t) = bc.forcing(t)

BoundaryStyle(::Type{<:AirTemperature}) = Dirichlet()

export AirTemperature
