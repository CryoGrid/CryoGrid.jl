abstract type GroundParameterization end

"""
    AbstractGround{Tpara<:GroundParameterization,Theat<:Optional{HeatBalance},Twater<:Optional{WaterBalance}} <: SubSurface

Base type for all ground layers defining heat and water balance schemes.
"""
abstract type AbstractGround{Tpara<:GroundParameterization,Theat<:Optional{HeatBalance},Twater<:Optional{WaterBalance}} <: SubSurface end

"""
    SoilParameterization

Base type for parameterizations of soil consituents.
"""
abstract type SoilParameterization <: GroundParameterization end

"""
    Soil{Tpara,Theat,Twater}

Type alias for any `AbstractGround` layer with a `SoilParameterization`.
"""
const Soil{Tpara,Theat,Twater} = AbstractGround{Tpara,Theat,Twater} where {Tpara<:SoilParameterization,Theat<:Optional{HeatBalance},Twater<:Optional{WaterBalance}}

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

"""
    Heterogeneous{Tpara,Taux} <: SoilParameterization

Special `SoilParameterization` which wraps another soil parameterization type
to indicate that it should be heterogeneous with over depth. Parameterizations
that support such configurations should provide dispatches for `Heterogeneous{...}`
that instantiate the relevant soil properties as on-grid state variables.
"""
Base.@kwdef struct Heterogeneous{Tpara,Taux} <: SoilParameterization
    para::Tpara
    aux::Taux = nothing
    Heterogeneous(para::SoilParameterization, aux=nothing) = new{typeof(para),typeof(aux)}(para, aux)
end

checksolver!(::SoilParameterization, ::Union{Nothing,SFCCSolver}) = nothing
checksolver!(::Heterogeneous, ::SFCCPreSolver) = error("SFCCPreSolver requires homogeneous soil properties in each layer.")
