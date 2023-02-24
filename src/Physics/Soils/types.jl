"""
    SoilParameterization

Base type for parameterizations of soil consituents.
"""
abstract type SoilParameterization end

"""
    Soil{Tpara,Tproc} <: SubSurface{Tproc}

Base type for all soil or soil-like layers.
"""
abstract type Soil{Tpara,Tproc} <: SubSurface{Tproc} end

"""
    HomogeneousSoil{Tpara<:SoilParameterization,Tprop,Tsp,TP} <: SubSurface{TP}

Generic, homogeneous Soil layer, i.e. material is assumed to be uniformly mixed.
"""
struct HomogeneousSoil{Tpara<:SoilParameterization,Tproc,Tsp,Tprop,Tsolver} <: Soil{Tpara,Tproc}
    para::Tpara # soil parameterization
    prop::Tprop # soil properties
    proc::Tproc # subsurface processes
    solver::Tsolver
    sp::Tsp # user-defined specialization
end
