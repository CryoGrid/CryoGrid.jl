const Dirichlet = CryoGrid.Dirichlet
const Neumann = CryoGrid.Neumann
# Boundary condition type aliases
const HeatBC = BoundaryProcess{T} where {HeatBalance<:T<:SubSurfaceProcess}
ConstantTemperature(value::UFloat"K") = ConstantBC(HeatBalance, Dirichlet, uconvert(u"°C", value))
ConstantTemperature(value) = ConstantBC(HeatBalance, Dirichlet, value)

# Boundary fluxes
function CryoGrid.boundaryflux(::Dirichlet, bc::HeatBC, top::Top, heat::HeatBalance, sub::SubSurface, stop, ssub)
    Δk = CryoGrid.thickness(sub, ssub, first) # using `thickness` allows for generic layer implementations
    @inbounds let Tupper=boundaryvalue(bc, stop),
        Tsub=ssub.T[1],
        k=ssub.k[1],
        δ=Δk/2; # distance to boundary
        Numerics.flux(Tupper, Tsub, δ, k)
    end
end
function CryoGrid.boundaryflux(::Dirichlet, bc::HeatBC, bot::Bottom, heat::HeatBalance, sub::SubSurface, sbot, ssub)
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
    TemperatureBC{E,F} <: BoundaryProcess{HeatBalance}

Represents a simple, Dirichlet temperature boundary condition for `HeatBalance` processes.
"""
struct TemperatureBC{E,F} <: BoundaryProcess{HeatBalance}
    T::F # boundary temperature
    effect::E # effect
    TemperatureBC(T::F, effect::E=nothing) where {F,E} = new{E,F}(T, effect)
end

CryoGrid.BCKind(::Type{<:TemperatureBC}) = Dirichlet()

CryoGrid.boundaryvalue(::TemperatureBC, state) = getscalar(state.T_ub)

CryoGrid.variables(::Top, ::TemperatureBC) = (
    Diagnostic(:T_ub, Scalar, u"K"),
)

CryoGrid.variables(::Bottom, ::TemperatureBC) = (
    Diagnostic(:T_lb, Scalar, u"K"),
)

function CryoGrid.computediagnostic!(::Top, bc::TemperatureBC, state)
    @setscalar state.T_ub = bc.T
end

function CryoGrid.computediagnostic!(::Bottom, bc::TemperatureBC, state)
    @setscalar state.T_lb = bc.T
end

Base.@kwdef struct NFactor{W,S} <: CryoGrid.BoundaryEffect
    nf::W = Param(1.0, domain=0..1) # applied when Tair <= 0
    nt::S = Param(1.0, domain=0..1) # applied when Tair > 0
end

CryoGrid.variables(::Top, bc::TemperatureBC{<:NFactor}) = (
    Diagnostic(:T_ub, Scalar, u"K"),
    Diagnostic(:nfactor, Scalar),
)

nfactor(Tair, nfw, nfs) = (Tair <= zero(Tair))*nfw + (Tair > zero(Tair))*nfs

function CryoGrid.computediagnostic!(::Top, bc::TemperatureBC{<:NFactor}, state)
    nfw = bc.effect.nf
    nfs = bc.effect.nt
    Tair = bc.T
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
    GroundHeatFlux(Qg::TQ, effect::TE=nothing) where {TQ,TE} = new{TE,TQ}(Qg, effect)
end

CryoGrid.boundaryvalue(bc::GroundHeatFlux, state) = bc.Qg

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

CryoGrid.BCKind(::Type{<:GeothermalHeatFlux}) = CryoGrid.Neumann()
