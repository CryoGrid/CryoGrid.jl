const Dirichlet = CryoGrid.Dirichlet
const Neumann = CryoGrid.Neumann
# Boundary condition type aliases
const HeatBC = BoundaryProcess{T} where {HeatBalance<:T<:SubSurfaceProcess}
ConstantTemperature(value::UFloat"K") = ConstantBC(HeatBalance, Dirichlet, uconvert(u"°C", value))
ConstantTemperature(value) = ConstantBC(HeatBalance, Dirichlet, value)

# Boundary fluxes
@inline function CryoGrid.boundaryflux(::Dirichlet, bc::HeatBC, top::Top, heat::HeatBalance, sub::SubSurface, stop, ssub)
    Δk = CryoGrid.thickness(sub, ssub, first) # using `thickness` allows for generic layer implementations
    @inbounds let Tupper=boundaryvalue(bc, stop),
        Tsub=ssub.T[1],
        k=ssub.k[1],
        δ=Δk/2; # distance to boundary
        Numerics.flux(Tupper, Tsub, δ, k)
    end
end
@inline function CryoGrid.boundaryflux(::Dirichlet, bc::HeatBC, bot::Bottom, heat::HeatBalance, sub::SubSurface, sbot, ssub)
    Δk = CryoGrid.thickness(sub, ssub, last) # using `thickness` allows for generic layer implementations
    @inbounds let Tlower=boundaryvalue(bc, sbot),
        Tsub=ssub.T[end],
        k=ssub.k[end],
        δ=Δk/2; # distance to boundary
        Numerics.flux(Tsub, Tlower, δ, k)
    end
end

function CryoGrid.interact!(top::Top, bc::HeatBC, sub::SubSurface, heat::HeatBalance, stop, ssub)
    ssub.jH[1] += boundaryflux(bc, top, heat, sub, stop, ssub)
    return nothing
end
function CryoGrid.interact!(sub::SubSurface, heat::HeatBalance, bot::Bottom, bc::HeatBC, ssub, sbot)
    ssub.jH[end] += boundaryflux(bc, bot, heat, sub, sbot, ssub)
    return nothing
end

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

CryoGrid.boundaryvalue(bc::TemperatureGradient, state) = getscalar(state.T_ub)

CryoGrid.variables(::Union{Top,Bottom}, bc::TemperatureGradient) = (
    Diagnostic(:T_ub, Scalar, u"K"),
)

function CryoGrid.updatestate!(::Top, bc::TemperatureGradient, state)
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

function CryoGrid.updatestate!(::Top, bc::TemperatureGradient{<:NFactor}, state)
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

CryoGrid.boundaryvalue(bc::GroundHeatFlux, state) = bc.Qg(state.t)

CryoGrid.BCKind(::Type{<:GroundHeatFlux}) = CryoGrid.Neumann()

"""
    GeothermalHeatFlux{TQ} <: BoundaryProcess{HeatBalance}

Represents a simple, forced Neumann heat flux boundary condition for `HeatBalance` processes.
"""
struct GeothermalHeatFlux{TQ} <: BoundaryProcess{HeatBalance}
    Qgeo::TQ
    GeothermalHeatFlux(Qgeo::TQ) where {TQ} = new{TQ}(Qgeo)
end

CryoGrid.boundaryvalue(bc::GeothermalHeatFlux, state) = -bc.Qgeo

CryoGrid.boundaryvalue(bc::GeothermalHeatFlux{<:Forcing}, state) = -bc.Qgeo(state.t)

CryoGrid.BCKind(::Type{<:GeothermalHeatFlux}) = CryoGrid.Neumann()
