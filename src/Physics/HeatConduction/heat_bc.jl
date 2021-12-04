# Boundary condition type aliases
const ConstantTemp = Constant{Dirichlet,Float"°C"}
ConstantTemp(value::UFloat"K") = Constant(Dirichlet, dustrip(u"°C", value))
ConstantTemp(value::UFloat"°C") = Constant(Dirichlet, dustrip(value))
const GeothermalHeatFlux = Constant{Neumann,Float"J/s/m^2"}
GeothermalHeatFlux(value::UFloat"J/s/m^2"=0.053xu"J/s/m^2") = Constant(Neumann, dustrip(value))

struct TemperatureGradient{E,F} <: BoundaryProcess
    T::F # temperature forcing
    effect::E # effect
    TemperatureGradient(T::F, effect::E=NoEffect()) where {F<:Forcing,E} = new{E,F}(T, effect)
end
BoundaryStyle(::Type{<:TemperatureGradient}) = Dirichlet()
BoundaryStyle(::Type{<:TemperatureGradient{<:Damping}}) = Neumann()

@with_kw struct NFactor{W,S} <: BoundaryEffect
    winterfactor::W = Param(1.0, bounds=(0.0,1.0)) # applied when Tair <= 0
    summerfactor::S = Param(1.0, bounds=(0.0,1.0)) # applied when Tair > 0
end

@inline (bc::TemperatureGradient)(l1,l2,p2,s1,s2) where {F} = bc.forcings.Tair(s1.t)
@inline function (bc::TemperatureGradient{<:NFactor})(l1,l2,p2::Heat,s1,s2)
    let nfw = bc.effect.winterfactor,
        nfs = bc.effect.summerfactor,
        Tair = bc.T(s1.t),
        nf = (Tair <= zero(Tair))*nfw + (Tair > zero(Tair))*nfs;
        nf*Tair # apply damping factor to air temperature
    end
end
@inline function (bc::TemperatureGradient{<:Damping})(l1,l2,p2::Heat,s1,s2)
    let Ttop = bc.T(s1.t),
        Tsub = s2.T[1],
        δ = Δ(s2.grids.k)[1], # grid cell size
        k_sub = s2.kc[1];
        bc.effect(Ttop, Tsub, k_sub, δ, s1.t)
    end
end
