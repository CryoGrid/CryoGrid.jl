struct TemperatureGradient{F} <: BoundaryProcess
    forcing::F
    TemperatureGradient(forcing::Forcing{Float"°C"}) = new{typeof(forcing)}(forcing)
end

@inline (bc::TemperatureGradient)(t) = bc.forcing(t)
@inline (bc::TemperatureGradient)(l1,l2,p2,s1,s2) = bc(s1.t)

BoundaryStyle(::Type{<:TemperatureGradient}) = Dirichlet()

struct NFactor{F,nf,nt} <: BoundaryProcess
    tgrad::TemperatureGradient{F}
    factor::Float64
    threshold::Float64
    NFactor(tgrad::TemperatureGradient{F}, factor::Float64=0.5, threshold::Float64=0.0, name::Symbol=:n) where {F} =
        new{F,Symbol(name,:_factor),Symbol(name,:_thresh)}(tgrad, factor, threshold)
end

@inline function (bc::NFactor{F,nf,nt})(l1,l2,p2,s1,s2) where {F,nf,nt}
    let factor = s1.params[nf] |> getscalar,
        thresh = s1.params[nt] |> getscalar,
        Tair = bc.tgrad(l1,l2,p2,s1,s2),
        check = Tair - thresh,
        factor_eff = (check <= zero(factor))*factor + (check > zero(factor));
        factor_eff*Tair
    end
end

variables(::Top, bc::NFactor{F,nf,nt}) where {F,nf,nt} = (Parameter(nf, bc.factor, 0..1), Parameter(nt, bc.threshold))
BoundaryStyle(::Type{<:NFactor}) = Dirichlet()

# Boundary condition type aliases
const ConstantTemp = Constant{Dirichlet,Float"°C"}
ConstantTemp(value::UFloat"K") = Constant{Dirichlet}(dustrip(u"°C", value))
ConstantTemp(value::UFloat"°C") = Constant{Dirichlet}(dustrip(value))
const GeothermalHeatFlux = Constant{Neumann,Float"J/s/m^2"}
GeothermalHeatFlux(value::UFloat"J/s/m^2"=0.053xu"J/s/m^2") = Constant{Neumann}(dustrip(value))
