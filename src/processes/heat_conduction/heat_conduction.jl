abstract type FreezeCurve end
struct FreeWater <: FreezeCurve end
# Default implementation of variables
variables(fc::FreezeCurve) = ()
export FreeWater, FreezeCurve

@with_kw struct HeatParams{T<:FreezeCurve} <: Params
    ρ::Float"kg/m^3" = 1000.0xu"kg/m^3" #[kg/m^3]
    Lsl::Float"J/kg" = 334000.0xu"J/kg" #[J/kg] (latent heat of fusion)
    L::Float"J/m^3" = (ρ*Lsl)xu"J/m^3" #[J/m^3] (specific latent heat of fusion)
    freezecurve::T = FreeWater()
end

"""
Alias and constructor for Profile specific to temperature.
"""
const TempProfile{D,Q,T} = Profile{D,1,Q,T} where {D,Q,T}
TempProfile(pairs::Pair{<:DistQuantity, <:TempQuantity}...) =
    Profile([d=>(uconvert(u"K",T),) for (d,T) in pairs]...;names=(:T,))

struct Heat{U,TParams} <: SubSurfaceProcess
    params::TParams
    profile::Union{Nothing,TempProfile}
    function Heat{stateunit}(profile::TProfile=nothing; kwargs...) where {stateunit, TProfile<:Union{Nothing,TempProfile}}
        @assert stateunit == u"J" || stateunit == u"K" "State unit type parameter must be either J or K"
        params = HeatParams(;kwargs...)
        new{stateunit,typeof(params)}(params,profile)
    end
end

Base.show(io::IO, h::Heat{U,P}) where {U,P} = print(io, "Heat{$U,$P}($(h.params))")

export Heat, HeatParams, TempProfile

freezecurve(heat::Heat) = heat.params.freezecurve
enthalpy(T::Real"K", C::Real"J/K/m^3", L::Real"J/m^3", θ::Real) = (T-273.15)*C + L*θ

export freezecurve, enthalpy

"""
    heatconduction!(T,ΔT,k,Δk,∂H)

1-D heat conduction/diffusion given T, k, and their deltas. Resulting enthalpy gradient is stored in ∂H.
Note that this function does not perform bounds checking. It is up to the user to ensure that all variables are
arrays of the correct length.
"""
function heatconduction!(T,ΔT,k,Δk,∂H)
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
"""
    boundaryflux(boundary::Boundary, bc::B, sub::SubSurface, h::Heat, sbound, ssub) where {B<:BoundaryProcess{Heat}}

Computes the flux dH/dt at the given boundary. Calls boundaryflux(BoundaryStyle(B),...) to allow for generic
implementations by boundary condition type.
"""
boundaryflux(boundary::Boundary, bc::B, sub::SubSurface, h::Heat, sbound, ssub) where {B<:BoundaryProcess{<:Heat}} =
    boundaryflux(BoundaryStyle(B), boundary, bc, sub, h, sbound, ssub)
boundaryflux(::Neumann, top::Top, bc::B, sub::SubSurface, h::Heat, stop, ssub) where {B<:BoundaryProcess{<:Heat}} =
    @inbounds let a = Δ(ssub.grids.k)[1]
        bc(top,sub,h,stop,ssub)/a
    end
boundaryflux(::Neumann, bot::Bottom, bc::B, sub::SubSurface, h::Heat, sbot, ssub) where {B<:BoundaryProcess{<:Heat}} =
    @inbounds let a = Δ(ssub.grids.k)[end]
        bc(bot,sub,h,sbot,ssub)/a
    end
function boundaryflux(::Dirichlet, top::Top, bc::B, sub::SubSurface, h::Heat, stop, ssub) where {B<:BoundaryProcess{<:Heat}}
    Δk = Δ(ssub.grids.k)
    @inbounds let Tupper=bc(top,sub,h,stop,ssub),
        Tsub=ssub.T[1],
        k=ssub.k[1],
        a=Δk[1],
        δ=(Δk[1]/2); # distance to surface
        -k*(Tsub-Tupper)/δ/a
    end
end
function boundaryflux(::Dirichlet, bot::Bottom, bc::B, sub::SubSurface, h::Heat, sbot, ssub) where {B<:BoundaryProcess{<:Heat}}
    Δk = Δ(ssub.grids.k)
    @inbounds let Tlower=bc(bot,sub,h,sbot,ssub),
        Tsub=ssub.T[end],
        k=ssub.k[end],
        a=Δk[end],
        δ=(Δk[end]/2); # distance to surface
        -k*(Tsub-Tlower)/δ/a
    end
end
"""
Generic top interaction. Computes flux dH at top cell.
"""
function interact!(top::Top, bc::B, sub::SubSurface, heat::Heat{u"J"}, stop, ssub) where {B<:BoundaryProcess{<:Heat}}
    @inbounds ssub.dH[1] += boundaryflux(top, bc, sub, heat, stop, ssub)
    return nothing # ensure no allocation
end
function interact!(top::Top, bc::B, sub::SubSurface, heat::Heat{u"K"}, stop, ssub) where {B<:BoundaryProcess{<:Heat}}
    @inbounds ssub.dH[1] += boundaryflux(top, bc, sub, heat, stop, ssub)
    @inbounds ssub.dT[1] += ssub.dH[1] / ssub.Ceff[1]
    return nothing # ensure no allocation
end
"""
Generic bottom interaction. Computes flux dH at bottom cell.
"""
function interact!(sub::SubSurface, heat::Heat{u"J"}, bot::Bottom, bc::B, ssub, sbot) where {B<:BoundaryProcess{<:Heat}}
    @inbounds ssub.dH[end] += boundaryflux(bot, bc, sub, heat, sbot, ssub)
    return nothing # ensure no allocation
end
function interact!(sub::SubSurface, heat::Heat{u"K"}, bot::Bottom, bc::B, ssub, sbot) where {B<:BoundaryProcess{<:Heat}}
    @inbounds ssub.dH[end] += boundaryflux(bot, bc, sub, heat, sbot, ssub)
    @inbounds ssub.dT[end] += ssub.dH[end] / ssub.Ceff[end]
    return nothing # ensure no allocation
end
# Free water freeze curve
"""
Implementation of "free water" freeze curve for any subsurface layer. Assumes that
'state' contains at least temperature (T), enthalpy (H), heat capacity (C),
total water content (θw), and liquid water content (θl).
"""
@inline function (fc::FreeWater)(layer::SubSurface, heat::Heat{u"J"}, state)
    @inline function enthalpyinv(H, C, L, θtot)
        let θtot = max(1.0e-8,θtot),
            Lθ = L*θtot,
            I_t = H > Lθ,
            I_f = H <= 0.0;
            T = (I_t*(H-Lθ) + I_f*H)/C + 273.15
        end
    end
    @inline function freezethaw(H, C, L, θtot)
        let θtot = max(1.0e-8,θtot),
            Lθ = L*θtot,
            I_t = H > Lθ,
            I_c = (H > 0.0) && (H <= Lθ);
            liquidfraction = I_c*(H/Lθ) + I_t
        end
    end
    L = heat.params.L
    @. state.T = enthalpyinv(state.H, state.C, L, state.θw)
    @. state.θl = freezethaw(state.H, state.C, L, state.θw)*state.θw
end
# Fallback (error) implementation for freeze curve
(fc::FreezeCurve)(layer::SubSurface, heat::Heat, state) =
    error("freeze curve $(typeof(fc)) not implemented for $(typeof(heat)) on layer $(typeof(layer))")

export heatconduction!, boundaryflux

include("soil/soilheat.jl")
