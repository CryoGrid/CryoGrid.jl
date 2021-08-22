# Boundary condition type aliases
const ConstantTemp = Constant{Dirichlet,Float"°C"}
ConstantTemp(value::UFloat"K") = Constant{Dirichlet}(dustrip(u"°C", value))
ConstantTemp(value::UFloat"°C") = Constant{Dirichlet}(dustrip(value))
const GeothermalHeatFlux = Constant{Neumann,Float"J/s/m^2"}
GeothermalHeatFlux(value::UFloat"J/s/m^2"=0.053xu"J/s/m^2") = Constant{Neumann}(dustrip(value))

abstract type TempGradTransform end
struct NoTransform <: TempGradTransform end
struct StationaryNFactor{nfw,nfs} <: TempGradTransform end
struct TwoStageNFactor{nfw1,nfw2,nfs1,nfs2,nftc} <: TempGradTransform end

struct TemperatureGradient{T,F} <: BoundaryProcess
    forcing::F
    transform::T
    TemperatureGradient(forcing::Forcing{Float"°C"}, transform::T=NoTransform()) where {T<:TempGradTransform} = new{T,typeof(forcing)}(forcing, transform)
end

NFactor(type::Symbol=:stationary, name::Symbol=:nf) = NFactor(Val{type}(), name)
NFactor(::Val{:stationary}, name::Symbol=:nf) = StationaryNFactor{Symbol(name,:w),Symbol(name,:s)}()
NFactor(::Val{:twostage}, name::Symbol=:nf) = TwoStageNFactor{Symbol(name,:w1),Symbol(name,:w2),Symbol(name,:s1),Symbol(name,:s2),Symbol(name,:tc)}()

@inline (bc::TemperatureGradient)(l1,l2,p2,s1,s2) where {F} = bc(s1.t)
@inline function (bc::TemperatureGradient{StationaryNFactor{nfw,nfs}})(l1,l2,p2,s1,s2) where {nfw,nfs}
    let nf_winter = s1.params[nfw] |> getscalar,
        nf_summer = s1.params[nfs] |> getscalar,
        Tair = bc.forcing(s1.t),
        factor_eff = (Tair <= zero(Tair))*nf_winter + (Tair > zero(Tair))*nf_summer;
        Tbc = factor_eff*Tair
    end
end
@inline function (bc::TemperatureGradient{TwoStageNFactor{nfw1,nfw2,nfs1,nfs2,nftc}})(l1,l2,p2,s1,s2) where {nfw1,nfw2,nfs1,nfs2,nftc}
    let nf_winter1 = s1.params[nfw1] |> getscalar,
        nf_winter2 = s1.params[nfw2] |> getscalar,
        nf_summer1 = s1.params[nfs1] |> getscalar,
        nf_summer2 = s1.params[nfs2] |> getscalar,
        nf_tchange = s1.params[nftc] |> getscalar,
        tcheck = s1.t >= nf_tchange,
        Tair = bc.forcing(s1.t),
        factor_eff = (Tair <= zero(Tair))*(1-tcheck)*nf_winter1 +
            (Tair > zero(Tair))*(1-tcheck)*nf_summer1 +
            (Tair <= zero(Tair))*tcheck*nf_winter2 +
            (Tair > zero(Tair))*tcheck*nf_summer2
        Tbc = factor_eff*Tair
    end
end

variables(::Top, bc::TemperatureGradient{StationaryNFactor{nfw,nfs}}) where {nfw,nfs} = (Parameter(nfw, 1.0, 0..1), Parameter(nfs, 1.0, 0..1))
variables(::Top, bc::TemperatureGradient{TwoStageNFactor{nfw1,nfw2,nfs1,nfs2,nftc}}) where {nfw1,nfw2,nfs1,nfs2,nftc} = (
    (Parameter(name, 1.0, 0..1) for name in (nfw1,nfw2,nfs1,nfs2))...,
    Parameter(nftc, 0.0, 0..Inf)
)

BoundaryStyle(::Type{<:TemperatureGradient}) = Dirichlet()
