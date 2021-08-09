struct TemperatureGradient{F} <: BoundaryProcess
    forcing::F
    TemperatureGradient(forcing::Forcing{Float"°C"}) = new{typeof(forcing)}(forcing)
end

@inline (bc::TemperatureGradient)(t) = bc.forcing(t)
@inline (bc::TemperatureGradient)(l1,l2,p2,s1,s2) = bc(s1.t)

BoundaryStyle(::Type{<:TemperatureGradient}) = Dirichlet()

struct NFactor{F,nf1,nf2,nt} <: BoundaryProcess
    tgrad::TemperatureGradient{F}
    factor1::Float64
    factor2::Float64
    threshold::Float64
    NFactor(tgrad::TemperatureGradient{F}, factor1::Float64=0.5, factor2::Float64=1.0, threshold::Float64=0.0, name::Symbol=:n) where {F} =
        new{F,Symbol(name,:_factor1),Symbol(name,:_factor2),Symbol(name,:_thresh)}(tgrad, factor1, factor2, threshold)
end

@inline function (bc::NFactor{F,nf1,nf2,nt})(l1,l2,p2,s1,s2) where {F,nf1,nf2,nt}
    let factor1 = s1.params[nf1] |> getscalar,
        factor2 = s1.params[nf2] |> getscalar,
        thresh = s1.params[nt] |> getscalar,
        Tair = bc.tgrad(l1,l2,p2,s1,s2),
        check = Tair - thresh,
        factor_eff = (check <= zero(check))*factor1 + (check > zero(check))*factor2;
        factor_eff*Tair
    end
end

variables(::Top, bc::NFactor{F,nf1,nf2,nt}) where {F,nf1,nf2,nt} = (Parameter(nf1, bc.factor1, 0..1), Parameter(nf2, bc.factor2, 0..1), Parameter(nt, bc.threshold))
BoundaryStyle(::Type{<:NFactor}) = Dirichlet()

# Boundary condition type aliases
const ConstantTemp = Constant{Dirichlet,Float"°C"}
ConstantTemp(value::UFloat"K") = Constant{Dirichlet}(dustrip(u"°C", value))
ConstantTemp(value::UFloat"°C") = Constant{Dirichlet}(dustrip(value))
const GeothermalHeatFlux = Constant{Neumann,Float"J/s/m^2"}
GeothermalHeatFlux(value::UFloat"J/s/m^2"=0.053xu"J/s/m^2") = Constant{Neumann}(dustrip(value))
