const Dirichlet = CryoGrid.Dirichlet
const Neumann = CryoGrid.Neumann
# Boundary condition type aliases
const HeatBC = BoundaryProcess{T} where {HeatBalance<:T<:SubSurfaceProcess}
ConstantTemperature(value::UFloat"K") = ConstantBC(HeatBalance, Dirichlet, uconvert(u"°C", value))
ConstantTemperature(value) = ConstantBC(HeatBalance, Dirichlet, value)
GeothermalHeatFlux(value=0.053u"W/m^2") = ConstantBC(HeatBalance, Neumann, value)

"""
    TemperatureGradient{E,F} <: BoundaryProcess{HeatBalance}

Represents a simple, forced Dirichlet temperature boundary condition for `HeatBalance` processes.
"""
struct TemperatureGradient{E,F} <: BoundaryProcess{HeatBalance}
    T::F # temperature forcing
    effect::E # effect
    TemperatureGradient(T::F, effect::E=nothing) where {F<:Forcing{u"°C"},E} = new{E,F}(T, effect)
end
CryoGrid.BCKind(::Type{<:TemperatureGradient}) = Dirichlet()
@inline CryoGrid.boundaryvalue(bc::TemperatureGradient, l1, ::HeatBalance, l2, s1, s2) = getscalar(s1.T_ub)

CryoGrid.variables(::Top, bc::TemperatureGradient) = (
    Diagnostic(:T_ub, Scalar, u"K"),
)
function CryoGrid.diagnosticstep!(::Top, bc::TemperatureGradient, state)
    @setscalar state.T_ub = bc.T(state.t)
end

Base.@kwdef struct NFactor{W,S} <: CryoGrid.BoundaryEffect
    nf::W = 1.0 # applied when Tair <= 0
    nt::S = 1.0 # applied when Tair > 0
end
CryoGrid.parameterize(nf::NFactor) = NFactor(
    nf = CryoGrid.parameterize(nf.nf, domain=0..1),
    nt = CryoGrid.parameterize(nf.nt, domain=0..1),
)
CryoGrid.variables(::Top, bc::TemperatureGradient{<:NFactor}) = (
    Diagnostic(:T_ub, Scalar, u"K"),
    Diagnostic(:nfactor, Scalar),
)
nfactor(Tair, nfw, nfs) = (Tair <= zero(Tair))*nfw + (Tair > zero(Tair))*nfs
function CryoGrid.diagnosticstep!(::Top, bc::TemperatureGradient{<:NFactor}, state)
    nfw = bc.effect.nf
    nfs = bc.effect.nt
    Tair = bc.T(state.t)
    @setscalar state.nfactor = nfactor(Tair, nfw, nfs)
    @setscalar state.T_ub = getscalar(state.nfactor)*Tair
end

"""
    GroundHeatFlux{TE,TQ} <: BoundaryProcess{HeatBalance}

Represents a simple, forced Neumann heat flux boundary condition for `HeatBalance` processes.
"""
struct GroundHeatFlux{TE,TQ} <: BoundaryProcess{HeatBalance}
	Qg::TQ
    effect::TE
    GroundHeatFlux(Qg::TQ, effect::TE=nothing) where {TQ<:Forcing{u"W/m^2"},TE} = new{TE,TQ}(Qg, effect)
end
CryoGrid.boundaryvalue(bc::GroundHeatFlux, ::Top, ::HeatBalance, ::SubSurface, stop, ssub) = bc.Qg(stop.t)
CryoGrid.BCKind(::Type{<:GroundHeatFlux}) = CryoGrid.Neumann()