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
    HomogeneousSoil{Tpara,Theat<:Optional{HeatBalance},Twater<:Optional{WaterBalance},Tsp} <: Soil{Tpara,Theat,Twater}

Generic representation of a soil layer.
"""
Base.@kwdef struct HomogeneousSoil{Tpara,Theat<:Optional{HeatBalance},Twater<:Optional{WaterBalance},Tsp} <: Soil{Tpara,Theat,Twater}
    para::Tpara = MineralOrganic() # soil parameterization
    heat::Theat = HeatBalance() # heat conduction
    water::Twater = nothing # water balance
    sp::Tsp = nothing # user-defined specialization
end
# Convenience constructor
HomogeneousSoil(para::SoilParameterization; kwargs...) = HomogeneousSoil(; para, kwargs...)
