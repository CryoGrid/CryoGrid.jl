# Boundary condition type aliases
const HeatBC = BoundaryProcess{T} where {Heat<:T<:SubSurfaceProcess}
const ConstantTemp = ConstantBC{Heat, Dirichlet,Float"°C"}
ConstantTemp(value::UFloat"K") = ConstantBC(Heat, Dirichlet, dustrip(u"°C", value))
ConstantTemp(value::UFloat"°C") = ConstantBC(Heat, Dirichlet, dustrip(value))
const GeothermalHeatFlux = ConstantBC{Heat, Neumann, Float"J/s/m^2"}
GeothermalHeatFlux(value::UFloat"J/s/m^2"=0.053xu"J/s/m^2") = ConstantBC(Heat, Neumann, dustrip(value))

struct TemperatureGradient{E,F} <: BoundaryProcess{Heat}
    T::F # temperature forcing
    effect::E # effect
    TemperatureGradient(T::F, effect::E=nothing) where {F<:Forcing,E} = new{E,F}(T, effect)
end
BoundaryStyle(::Type{<:TemperatureGradient}) = Dirichlet()
@inline boundaryvalue(bc::TemperatureGradient, l1, ::Heat, l2, s1, s2) = getscalar(s1.T_ub)

variables(::Top, bc::TemperatureGradient) = (
    Diagnostic(:T_ub, Scalar, u"K"),
)
function diagnosticstep!(::Top, bc::TemperatureGradient, state)
    @setscalar state.T_ub = bc.T(state.t)
end

Base.@kwdef struct NFactor{W,S} <: BoundaryEffect
    nf::W = Param(1.0, bounds=(0.0,1.0)) # applied when Tair <= 0
    nt::S = Param(1.0, bounds=(0.0,1.0)) # applied when Tair > 0
end
variables(::Top, bc::TemperatureGradient{<:NFactor}) = (
    Diagnostic(:T_ub, Scalar, u"K"),
    Diagnostic(:nfactor, Scalar),
)
nfactor(Tair, nfw, nfs) = (Tair <= zero(Tair))*nfw + (Tair > zero(Tair))*nfs
function diagnosticstep!(::Top, bc::TemperatureGradient{<:NFactor}, state)
    nfw = bc.effect.nf
    nfs = bc.effect.nt
    Tair = bc.T(state.t)
    @setscalar state.nfactor = nfactor(Tair, nfw, nfs)
    @setscalar state.T_ub = getscalar(state.nfactor)*Tair
end
