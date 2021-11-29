# Boundary condition type aliases
const ConstantTemp = Constant{Dirichlet,Float"°C"}
ConstantTemp(value::UFloat"K") = Constant(Dirichlet, dustrip(u"°C", value))
ConstantTemp(value::UFloat"°C") = Constant(Dirichlet, dustrip(value))
const GeothermalHeatFlux = Constant{Neumann,Float"J/s/m^2"}
GeothermalHeatFlux(value::UFloat"J/s/m^2"=0.053xu"J/s/m^2") = Constant(Neumann, dustrip(value))

abstract type TemperatureGradientEffect end
struct NoEffect <: TemperatureGradientEffect end
@with_kw struct NFactor{W,S} <: TemperatureGradientEffect
    winterfactor::W = Param(1.0, bounds=(0.0,1.0)) # applied when Tair <= 0
    summerfactor::S = Param(1.0, bounds=(0.0,1.0)) # applied when Tair > 0
end
@with_kw struct SnowDamping{P,R} <: TemperatureGradientEffect
    p::P = Param(0.0, bounds=(0.0,1.0))
    r::R = Param(1.0, bounds=(0.0,1.0))
end

struct TemperatureGradient{E,F} <: BoundaryProcess
    forcings::F
    effect::E
    TemperatureGradient(forcings::F, effect::E=NoEffect()) where {F<:Forcings,E<:TemperatureGradientEffect} = new{E,F}(forcings, effect)
    function TemperatureGradient(Tair::Forcing, effect::E=NoEffect()) where {E<:TemperatureGradientEffect}
        forcings = Forcings(Tair=Tair)
        new{E,typeof(forcings)}(forcings, effect)
    end
end
BoundaryStyle(::Type{<:TemperatureGradient}) = Dirichlet()

@inline (bc::TemperatureGradient)(l1,l2,p2,s1,s2) where {F} = bc.forcings.Tair(s1.t)
@inline function (bc::TemperatureGradient{<:NFactor})(l1,l2,p2,s1,s2)
    let nfw = bc.effect.winterfactor,
        nfs = bc.effect.summerfactor,
        Tair = bc.forcings.Tair(s1.t),
        nf = (Tair <= zero(Tair))*nfw + (Tair > zero(Tair))*nfs;
        nf*Tair # apply damping factor to air temperature
    end
end
@inline function (bc::TemperatureGradient{<:SnowDamping})(l1,l2,p2,s1,s2)
    let p = bc.effect.p,
        r = bc.effect.r,
        sd = bc.forcings.Dsn(s1.t),
        Tair = bc.forcings.Tair(s1.t),
        nf = p + (1-p)exp(-sd/r);
        nf*Tair # apply damping factor to air temperature
    end
end
