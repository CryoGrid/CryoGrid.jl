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
    function SimpleSoil(para, heat, water, aux)
        solver = Heat.fcsolver(heat)
        @assert checksolver(para, solver) "SFCCPreSolver requires homogeneous soil properties in each layer."
        new{typeof(para),typeof(heat),typeof(water),typeof(aux)}(para, heat, water, aux)
    end
end
# Convenience constructor
SimpleSoil(para::SoilParameterization; kwargs...) = SimpleSoil(; para, kwargs...)

"""
    Heterogeneous{Tpara,Taux} <: SoilParameterization

Specialized `SoilParameterization` which wraps another soil parameterization
to make it heterogeneous with respect to depth. Parameterizations which support
this should provide dispatches for `Heterogeneous{...}` that instantiate the
relevant soil properties as on-grid state variables.
"""
Base.@kwdef struct Heterogeneous{Tpara,Taux} <: SoilParameterization
    para::Tpara
    aux::Taux = nothing
    Heterogeneous(para::SoilParameterization, aux=nothing) = new{typeof(para),typeof(aux)}(para, aux)
end

checksolver(::SoilParameterization, ::Union{Nothing,SFCCSolver}) = true
checksolver(::Heterogeneous, ::SFCCPreSolver) = false
