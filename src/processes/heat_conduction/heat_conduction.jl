abstract type FreezeCurve end
struct LinearFC <: FreezeCurve end

@with_kw struct HeatParams{T} <: Params @deftype
    ρ::UFloat"kg/m^3" = 1000.0U"kg/m^3" #[kg/m^3] (default value assumes pure water)
    Lsl::UFloat"J/kg" = 334000.0U"J/kg" #[J/kg] (latent heat of fusion)
end

HeatParams(::Type{TFreezeCurve}=LinearFC;kwargs...) where {TFreezeCurve} = HeatParams{TFreezeCurve}(kwargs...)

const TempProfile{D,Q} = Profile{D,1,Q} where {D,Q<:DistQuantity}
TempProfile(pairs::Pair{UFloat"m", UFloat"°C"}...) = Profile([d=>(T,) for (d,T) in pairs]...;names=(:T,))

struct Heat{U,TParams,TProfile} <: SubSurfaceProcess
    params::TParams
    profile::TProfile
    Heat{UType"J"}(profile::TProfile, params::HeatParams=HeatParams()) where {TProfile<:TempProfile} =
        new{UType"J",typeof(params),TProfile}(params)
    Heat{UType"K"}(profile::TProfile, params::HeatParams=HeatParams()) where {TProfile<:TempProfile} =
        new{UType"K",typeof(params),TProfile}(params)
end

include("hc_soil.jl")

export Heat, HeatConfig, FreezeCurve, LinearFC
export enthalpy, enthalpyInv
