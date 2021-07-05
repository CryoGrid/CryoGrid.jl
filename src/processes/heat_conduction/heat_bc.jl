struct TemperatureGradient{F} <: BoundaryProcess{Heat}
    forcing::F
    TemperatureGradient(forcing::Forcing{Float"°C"}) = new{typeof(forcing)}(forcing)
end

@inline (bc::TemperatureGradient)(t) = bc.forcing(t)
@inline (bc::TemperatureGradient)(l1,l2,p2,s1,s2) = bc(s1.t)

CryoGrid.BoundaryStyle(::Type{<:TemperatureGradient}) = Dirichlet()

export TemperatureGradient

# Boundary condition type aliases
const ConstantTemp = Constant{Heat,Dirichlet,Float"K"}
ConstantTemp(value::UFloat"K") = Constant{Heat,Dirichlet}(dustrip(value))
ConstantTemp(value::UFloat"°C") = Constant{Heat,Dirichlet}(dustrip(u"K",value))
const GeothermalHeatFlux = Constant{Heat,Neumann,Float"J/s/m^2"}
GeothermalHeatFlux(value::UFloat"J/s/m^2"=0.053xu"J/s/m^2") = Constant{Heat,Neumann}(dustrip(value))

export ConstantTemp, GeothermalHeatFlux
