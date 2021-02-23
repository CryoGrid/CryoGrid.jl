abstract type FreezeCurve end
struct FreeWater <: FreezeCurve end

export FreezeCurve, FreeWater

@with_kw struct HeatParams{T} <: Params
    ρ::Float"kg/m^3" = 1000.0xu"kg/m^3" #[kg/m^3] (default value assumes pure water)
    Lsl::Float"J/kg" = 334000.0xu"J/kg" #[J/kg] (latent heat of fusion)
end

"""
Alias and constructor for Profile specific to temperature.
"""
const TempProfile{D,Q,T} = Profile{D,1,Q,T} where {D,Q,T}
TempProfile(pairs::Pair{<:DistQuantity, <:TempQuantity}...) =
    Profile([d=>(uconvert(u"K",T),) for (d,T) in pairs]...;names=(:T,))

struct Heat{U,TParams,TProfile} <: SubSurfaceProcess
    params::TParams
    profile::TProfile
    Heat{UT"J"}(profile::TProfile, params::HeatParams=HeatParams{FreeWater}()) where {TProfile<:TempProfile} =
        new{UT"J",typeof(params),TProfile}(params,profile)
    Heat{UT"K"}(profile::TProfile, params::HeatParams=HeatParams{FreeWater}()) where {TProfile<:TempProfile} =
        new{UT"K",typeof(params),TProfile}(params,profile)
end

ρ(heat::Heat) = heat.params.ρ
Lsl(heat::Heat) = heat.params.Lsl

export Heat, HeatParams, TempProfile, ρ, Lsl

# Boundary condition type aliases
const ConstantAirTemp = Constant{Heat,Dirichlet,Float"K"}
ConstantAirTemp(value::UFloat"K") = Constant{Heat,Dirichlet}(dustrip(value))
ConstantAirTemp(value::UFloat"°C") = Constant{Heat,Dirichlet}(dustrip(u"K",value))
const GeothermalHeatFlux = Constant{Heat,Neumann,Float"J/s"}
GeothermalHeatFlux(value::UFloat"J/s") = Constant{Heat,Neumann}(dustrip(value))

export ConstantAirTemp, GeothermalHeatFlux

include("hc_soil_H.jl")
