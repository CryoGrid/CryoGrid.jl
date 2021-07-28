struct TemperatureGradient{F} <: BoundaryProcess
    forcing::F
    TemperatureGradient(forcing::Forcing{Float"°C"}) = new{typeof(forcing)}(forcing)
end

@inline (bc::TemperatureGradient)(t) = bc.forcing(t)
@inline (bc::TemperatureGradient)(l1,l2,p2,s1,s2) = bc(s1.t)

BoundaryStyle(::Type{<:TemperatureGradient}) = Dirichlet()

struct NFactor{nf,nt} <: BoundaryProcess
    factor::Float64
    threshold::Float64
    NFactor(factor::Float64=0.5, threshold::Float64=0.0, name::Symbol=:n) = new{Symbol(name,:_factor),Symbol(name,:_thresh)}(factor, threshold)
end

@inline function (bc::NFactor{nf,nt})(l1,l2,p2,s1,s2) where {nf,nt}
    let factor = s1.params[nf] |> getscalar,
        thresh = s1.params[nt] |> getscalar,
        check = s2.T[1] - thresh;
        # if condition passes: 1*factor + 0
        # otherwise: 0*factor + 1
        (check <= zero(factor))*factor + (check > zero(factor))
    end
end

variables(::Top, bc::NFactor{nf,nt}) where {nf,nt} = (Parameter(nf, bc.factor, 0..1), Parameter(nt, bc.threshold))
BoundaryStyle(::Type{<:NFactor}) = Dirichlet()

# Boundary condition type aliases
const ConstantTemp = Constant{Dirichlet,Float"°C"}
ConstantTemp(value::UFloat"K") = Constant{Dirichlet}(dustrip(u"°C", value))
ConstantTemp(value::UFloat"°C") = Constant{Dirichlet}(dustrip(value))
const GeothermalHeatFlux = Constant{Neumann,Float"J/s/m^2"}
GeothermalHeatFlux(value::UFloat"J/s/m^2"=0.053xu"J/s/m^2") = Constant{Neumann}(dustrip(value))
