abstract type GroundParameterization end

"""
    AbstractGround{Tpara<:GroundParameterization,Theat<:HeatBalance,Twater<:WaterBalance} <: SubSurface

Base type for all ground layers defining heat and water balance schemes.
"""
abstract type AbstractGround{Tpara<:GroundParameterization,Theat<:Optional{HeatBalance},Twater<:Optional{WaterBalance}} <: SubSurface end

"""
    Ground{Tpara,Theat<:HeatBalance,Twater<:WaterBalance,Taux} <: Soil{Tpara,Theat,Twater}

Generic representation of a `Ground` layer with material parameterization `para`.
"""
Base.@kwdef struct Ground{Tpara,Theat<:Optional{HeatBalance},Twater<:Optional{WaterBalance},Tsolver,Taux} <: AbstractGround{Tpara,Theat,Twater}
    para::Tpara = SimpleSoil() # ground parameterization
    heat::Theat = HeatBalance() # heat conduction
    water::Twater = nothing # water balance
    fcsolver::Tsolver = default_fcsolver(para, heat, water)
    aux::Taux = nothing # user-defined specialization
end

# Convenience constructors
Ground(para::GroundParameterization; kwargs...) = Ground(; para, kwargs...)

default_fcsolver(::GroundParameterization, ::Nothing, ::WaterBalance) = nothing
default_fcsolver(::GroundParameterization, ::HeatBalance, ::Nothing) = SFCCPreSolver(FreezeCurves.SFCCPreSolverCache1D())
default_fcsolver(::GroundParameterization, ::HeatBalance, ::WaterBalance) = SFCCPreSolver(FreezeCurves.SFCCPreSolverCacheND())

fcsolver(ground::Ground) = ground.fcsolver

# CryoGrid methods

CryoGrid.processes(g::Ground{<:GroundParameterization,<:HeatBalance,Nothing}) = g.heat
CryoGrid.processes(g::Ground{<:GroundParameterization,Nothing,<:WaterBalance}) = g.water
CryoGrid.processes(g::Ground{<:GroundParameterization,<:HeatBalance,<:WaterBalance}) = Coupled(g.water, g.heat)
