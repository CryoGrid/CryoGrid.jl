"""
    SoilParameterization

Base type for parameterizations of soil consituents.
"""
abstract type SoilParameterization end

"""
    Soil{Tpara<:SoilParameterization,Theat<:Optional{HeatBalance},Twater<:Optional{WaterBalance}} <: SubSurface

Base type for all soil or soil-like layers.
"""
abstract type Soil{Tpara<:SoilParameterization,Theat<:Optional{HeatBalance},Twater<:Optional{WaterBalance}} <: SubSurface end

"""
    SimpleSoil{Tpara,Theat<:Optional{HeatBalance},Twater<:Optional{WaterBalance},Taux} <: Soil{Tpara,Theat,Twater}

Generic representation of a soil layer.
"""
Base.@kwdef struct SimpleSoil{Tpara,Theat<:Optional{HeatBalance},Twater<:Optional{WaterBalance},Taux} <: Soil{Tpara,Theat,Twater}
    para::Tpara = MineralOrganic() # soil parameterization
    heat::Theat = HeatBalance() # heat conduction
    water::Twater = nothing # water balance
    aux::Taux = nothing # user-defined specialization
end
# Convenience constructor
SimpleSoil(para::SoilParameterization; kwargs...) = SimpleSoil(; para, kwargs...)
