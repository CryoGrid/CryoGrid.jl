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
Base.@kwdef struct Ground{Tpara,Theat<:Optional{HeatBalance},Twater<:Optional{WaterBalance},Taux} <: AbstractGround{Tpara,Theat,Twater}
    para::Tpara = MineralOrganic() # ground parameterization
    heat::Theat = HeatBalance() # heat conduction
    water::Twater = nothing # water balance
    aux::Taux = nothing # user-defined specialization
    function Ground(para, heat, water, aux)
        solver = Heat.fcsolver(heat)
        checksolver!(para, solver) 
        new{typeof(para),typeof(heat),typeof(water),typeof(aux)}(para, heat, water, aux)
    end
end
# Convenience constructors
Ground(para::GroundParameterization; kwargs...) = Ground(; para, kwargs...)
Ground(knot::ProfileKnot{T,<:GroundParameterization}; kwargs...) where {T} = Ground(; para=knot.value, kwargs...)
