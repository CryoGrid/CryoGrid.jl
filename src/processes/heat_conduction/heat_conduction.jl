abstract type FreezeCurve end
struct LinearFC <: FreezeCurve end

export FreezeCurve, LinearFC

@with_kw struct HeatParams{T} <: Params
    ρ::UFloat"kg/m^3" = 1000.0U"kg/m^3" #[kg/m^3] (default value assumes pure water)
    Lsl::UFloat"J/kg" = 334000.0U"J/kg" #[J/kg] (latent heat of fusion)
end

#HeatParams(::Type{TFreezeCurve}=LinearFC;kwargs...) where {TFreezeCurve} = HeatParams{TFreezeCurve}(kwargs...)

const TempProfile{D,Q,T} = Profile{D,1,Q,T} where {D,Q,T}
TempProfile(pairs::Pair{<:DistQuantity, <:TempQuantity}...) =
    Profile([d=>(uconvert(u"K",T),) for (d,T) in pairs]...;names=(:T,))

struct Heat{U,TParams,TProfile} <: SubSurfaceProcess
    params::TParams
    profile::TProfile
    Heat{UT"J"}(profile::TProfile, params::HeatParams=HeatParams{LinearFC}()) where {TProfile<:TempProfile} =
        new{UT"J",typeof(params),TProfile}(params,profile)
    Heat{UT"K"}(profile::TProfile, params::HeatParams=HeatParams{LinearFC}()) where {TProfile<:TempProfile} =
        new{UT"K",typeof(params),TProfile}(params,profile)
end

ρ(heat::Heat) = heat.params.ρ
Lsl(heat::Heat) = heat.params.Lsl

export Heat, HeatParams, TempProfile, ρ, Lsl

# Boundary condition type aliases
const ConstantAirTemp = Constant{Heat,Dirichlet,UFloat"K"}
ConstantAirTemp(value::UFloat"K") = Constant{Heat,Dirichlet}(value)
const GeothermalHeatFlux = Constant{Heat,Neumann,UFloat"J"}
GeothermalHeatFlux(value::UFloat"J") = Constant{Heat,Neumann}(value)

export ConstantAirTemp, GeothermalHeatFlux

include("hc_soil.jl")
