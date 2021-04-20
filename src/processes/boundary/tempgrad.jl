struct TemperatureGradient{F} <: BoundaryProcess{Heat}
    forcing::F
    TemperatureGradient(forcing::Forcing{Float"°C"}) = new{typeof(forcing)}(forcing)
end

@inline (bc::TemperatureGradient)(t) = bc.forcing(t)
@inline (bc::TemperatureGradient)(l1,l2,p2,s1,s2) = bc(s1.t)

BoundaryStyle(::Type{<:TemperatureGradient}) = Dirichlet()

export TemperatureGradient
