struct Constant{P,S,T} <: BoundaryProcess{P}
    value::T
    Constant{P,S}(value::T) where {P<:SubSurfaceProcess,S<:BoundaryStyle,T} = new{P,S,T}(value)
end

# Arguments are irrelevant for Constant, so we can just use args...
(bc::Constant)(args...) = bc.value

BoundaryStyle(::Type{<:Constant{P,S}}) where {P,S} = S()

export Constant

# Boundary condition type aliases
const ConstantAirTemp = Constant{Heat,Dirichlet,Float"K"}
ConstantAirTemp(value::UFloat"K") = Constant{Heat,Dirichlet}(dustrip(value))
ConstantAirTemp(value::UFloat"°C") = Constant{Heat,Dirichlet}(dustrip(u"K",value))
const GeothermalHeatFlux = Constant{Heat,Neumann,Float"J/s"}
GeothermalHeatFlux(value::UFloat"J/s") = Constant{Heat,Neumann}(dustrip(value))

export ConstantAirTemp, GeothermalHeatFlux
