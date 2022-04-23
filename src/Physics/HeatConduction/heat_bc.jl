# Boundary condition type aliases
const HeatBC = BoundaryProcess{T} where {Heat<:T<:BoundaryProcess}
const ConstantTemp = ConstantBC{Heat, Dirichlet,Float"°C"}
ConstantTemp(value::UFloat"K") = ConstantBC(Heat, Dirichlet, dustrip(u"°C", value))
ConstantTemp(value::UFloat"°C") = ConstantBC(Heat, Dirichlet, dustrip(value))
const GeothermalHeatFlux = ConstantBC{Heat, Neumann,Float"J/s/m^2"}
GeothermalHeatFlux(value::UFloat"J/s/m^2"=0.053xu"J/s/m^2") = ConstantBC(Heat, Neumann, dustrip(value))

struct TemperatureGradient{E,F} <: BoundaryProcess{Heat}
    T::F # temperature forcing
    effect::E # effect
    TemperatureGradient(T::F, effect::E=nothing) where {F<:Forcing,E} = new{E,F}(T, effect)
end
BoundaryStyle(::Type{<:TemperatureGradient}) = Dirichlet()
BoundaryStyle(::Type{<:TemperatureGradient{<:Damping}}) = Neumann()
@inline boundaryvalue(bc::TemperatureGradient, l1, p2, l2, s1, s2) = bc.T(s1.t)

Base.@kwdef struct NFactor{W,S} <: BoundaryEffect
    nf::W = Param(1.0, bounds=(0.0,1.0)) # applied when Tair <= 0
    nt::S = Param(1.0, bounds=(0.0,1.0)) # applied when Tair > 0
end
@inline function boundaryvalue(bc::TemperatureGradient{<:NFactor}, l1, ::Heat, l2, s1, s2)
    let nfw = bc.effect.nf,
        nfs = bc.effect.nt,
        Tair = bc.T(s1.t),
        nf = (Tair <= zero(Tair))*nfw + (Tair > zero(Tair))*nfs;
        nf*Tair # apply damping factor to air temperature
    end
end
# damped temperature gradient
@inline function boundaryvalue(bc::TemperatureGradient{<:Damping}, l1, ::Heat, l2, s1,s2)
    let Ttop = bc.T(s1.t),
        Tsub = s2.T[1],
        δ = Δ(s2.grids.k)[1], # grid cell size
        k_sub = s2.kc[1];
        bc.effect(Ttop, Tsub, k_sub, δ, s1.t)
    end
end
