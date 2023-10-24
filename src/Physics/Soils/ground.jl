abstract type GroundParameterization end

"""
    AbstractGround{Tpara<:GroundParameterization,Theat<:Optional{HeatBalance},Twater<:Optional{WaterBalance}} <: SubSurface

Base type for all ground layers defining heat and water balance schemes.
"""
abstract type AbstractGround{Tpara<:GroundParameterization,Theat<:Optional{HeatBalance},Twater<:Optional{WaterBalance}} <: SubSurface end

"""
    Ground{Tpara,Theat<:Optional{HeatBalance},Twater<:Optional{WaterBalance},Taux} <: Soil{Tpara,Theat,Twater}

Generic representation of a `Ground` layer with material parameterization `para`.
"""
Base.@kwdef struct Ground{Tpara,Theat<:Optional{HeatBalance},Twater<:Optional{WaterBalance},Tsolver,Taux} <: AbstractGround{Tpara,Theat,Twater}
    para::Tpara = MineralOrganic() # ground parameterization
    heat::Theat = HeatBalance() # heat conduction
    water::Twater = nothing # water balance
    fcsolver::Tsolver = default_fcsolver(heat, water)
    aux::Taux = nothing # user-defined specialization
end
# Convenience constructors
Ground(para::GroundParameterization; kwargs...) = Ground(; para, kwargs...)
Ground(knot::ProfileKnot{T,<:GroundParameterization}; kwargs...) where {T} = Ground(; para=knot.value, kwargs...)

default_fcsolver(::Any, ::Any) = nothing
default_fcsolver(::HeatBalance{<:SFCC}, ::Nothing) = SFCCPreSolver(FreezeCurves.SFCCPreSolverCache1D())
default_fcsolver(::HeatBalance{<:SFCC}, ::WaterBalance) = SFCCPreSolver(FreezeCurves.SFCCPreSolverCacheND())

fcsolver(ground::Ground) = ground.fcsolver
