struct Constant{P,S,T} <: BoundaryProcess{P}
    value::T
    Constant{P,S}(value::T) where {P<:SubSurfaceProcess,S<:BoundaryStyle,T} = new{P,S,T}(value)
end

# Arguments are irrelevant for Constant, so we can just use args...
(bc::Constant)(args...) = bc.value

BoundaryStyle(::Type{<:Constant{P,S}}) where {P,S} = S()

export Constant

# Boundary condition type aliases
const ConstantTemp = Constant{Heat,Dirichlet,Float"K"}
ConstantTemp(value::UFloat"K") = Constant{Heat,Dirichlet}(dustrip(value))
ConstantTemp(value::UFloat"°C") = Constant{Heat,Dirichlet}(dustrip(u"K",value))
const GeothermalHeatFlux = Constant{Heat,Neumann,Float"J/s/m^2"}
GeothermalHeatFlux(value::UFloat"J/s/m^2"=0.053xu"J/s/m^2") = Constant{Heat,Neumann}(dustrip(value))

export ConstantTemp, GeothermalHeatFlux
