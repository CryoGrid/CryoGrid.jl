"""
    Heterogeneous{V,N,D,Taux} <: SoilParameterization

Special `SoilParameterization` which wraps a `Profile` of another soil parameterization type
to indicate that it should be heterogeneous with over depth.
"""
Base.@kwdef struct Heterogeneous{V,N,IT,VT,Taux} <: SoilParameterization
    profile::SoilProfile{N,IT,VT}
    aux::Taux = nothing
    Heterogeneous(::SoilProfile) = error("SoilProfile for heterogeneous layer must have uniform parameterization types (but the parameters may vary).")
    Heterogeneous(profile::SoilProfile{N,IT,VT}, aux=nothing) where {N,V<:SoilParameterization,IT<:NTuple{N},VT<:NTuple{N,V}} = new{V,N,IT,VT,typeof(aux)}(profile, aux)
end

"""
    Ground(soilprofile::SoilProfile; kwargs...)
"""
Ground(soilprofile::SoilProfile; kwargs...) = Ground(Heterogeneous(soilprofile); kwargs...)

# add dispatch for default_fcsolver that selects the ND presolver
default_fcsolver(::Heterogeneous, ::HeatBalance, ::Any) = SFCCPreSolver(FreezeCurves.SFCCPreSolverCacheND())

include("simple.jl")

include("surfex.jl")
