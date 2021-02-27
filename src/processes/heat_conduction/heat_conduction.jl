abstract type PhaseChangeStyle end
struct InstantFC <: PhaseChangeStyle end
struct FreeWaterFC <: PhaseChangeStyle end

export PhaseChangeStyle, FreeWaterFC, InstantFC, NoPhaseChange

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
    Heat{UT"J"}(profile::TProfile, params::HeatParams=HeatParams{FreeWaterFC}()) where {TProfile<:TempProfile} =
        new{UT"J",typeof(params),TProfile}(params,profile)
    Heat{UT"K"}(profile::TProfile, params::HeatParams=HeatParams{FreeWaterFC}()) where {TProfile<:TempProfile} =
        new{UT"K",typeof(params),TProfile}(params,profile)
end

ρ(heat::Heat) = heat.params.ρ
Lsl(heat::Heat) = heat.params.Lsl
"""
    heatconduction(T,ΔT,k,Δk,∂H)

1-D heat conduction/diffusion given T, k, and their deltas. Resulting enthalpy gradient is stored in ∂H.
Note that this function does not perform bounds checking. It is up to the user to ensure that all variables are
arrays of the correct length.
"""
function heatconduction(T,ΔT,k,Δk,∂H)
    # upper boundary
    @inbounds ∂H[1] += let T₂=T[2],
        T₁=T[1],
        k=k[2],
        δ=ΔT[1],
        a=Δk[1];
        k*(T₂-T₁)/δ/a
    end
    # diffusion on non-boundary cells
    @inbounds let T = T,
        k = (@view k[2:end-1]),
        Δk = (@view Δk[2:end-1]),
        ∂H = (@view ∂H[2:end-1]);
        ∇²(T, ΔT, k, Δk, ∂H)
    end
    # lower boundary
    @inbounds ∂H[end] += let T₂=T[end],
        T₁=T[end-1],
        k=k[end-1],
        δ=ΔT[end],
        a=Δk[end];
        -k*(T₂-T₁)/δ/a
    end
    return nothing
end

export Heat, HeatParams, TempProfile, ρ, Lsl, heatconduction

# Boundary condition type aliases
const ConstantAirTemp = Constant{Heat,Dirichlet,Float"K"}
ConstantAirTemp(value::UFloat"K") = Constant{Heat,Dirichlet}(dustrip(value))
ConstantAirTemp(value::UFloat"°C") = Constant{Heat,Dirichlet}(dustrip(u"K",value))
const GeothermalHeatFlux = Constant{Heat,Neumann,Float"J/s"}
GeothermalHeatFlux(value::UFloat"J/s") = Constant{Heat,Neumann}(dustrip(value))

export ConstantAirTemp, GeothermalHeatFlux

include("airtemp.jl")
include("hc_soil_H.jl")
