# Boundary condition type aliases
const ConstantTemp = Constant{Dirichlet,Float"°C"}
ConstantTemp(value::UFloat"K") = Constant(Dirichlet, dustrip(u"°C", value))
ConstantTemp(value::UFloat"°C") = Constant(Dirichlet, dustrip(value))
const GeothermalHeatFlux = Constant{Neumann,Float"J/s/m^2"}
GeothermalHeatFlux(value::UFloat"J/s/m^2"=0.053xu"J/s/m^2") = Constant(Neumann, dustrip(value))

abstract type TempGradTransform end
struct NoTransform <: TempGradTransform end
@with_kw struct NFactor{T} <: TempGradTransform
    winterfactor::T = Param(1.0, bounds=(0.0,1.0))
    summerfactor::T = Param(1.0, bounds=(0.0,1.0))
end

struct TemperatureGradient{T,F} <: BoundaryProcess
    forcings::F
    transform::T
    TemperatureGradient(forcings::F, transform::T=NoTransform()) where {F<:Forcings,T<:TempGradTransform} = new{T,F}(forcings, transform)
    function TemperatureGradient(Tair::Forcing, transform::T=NoTransform()) where {T<:TempGradTransform}
        forcings = Forcings(Tair=Tair)
        new{T,typeof(forcings)}(forcings, transform)
    end
end
BoundaryStyle(::Type{<:TemperatureGradient}) = Dirichlet()

@inline (bc::TemperatureGradient)(l1,l2,p2,s1,s2) where {F} = bc.forcings.Tair(s1.t)
@inline function (bc::TemperatureGradient{<:NFactor})(l1,l2,p2,s1,s2)
    let nfw = bc.transform.winterfactor,
        nfs = bc.transform.summerfactor,
        Tair = bc.forcings.Tair(s1.t),
        factor_eff = (Tair <= zero(Tair))*nfw + (Tair > zero(Tair))*nfs;
        Tbc = factor_eff*Tair
    end
end
