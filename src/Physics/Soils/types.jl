"""
    SoilParameterization

Base type for parameterizations of soil consituents.
"""
abstract type SoilParameterization end

"""
Generic container for numerical constants related to soil processes.
"""
Utils.@properties SoilProperties()

"""
    SoilProperties(para::SoilParameterization, proc::Process)

Constructor for `SoilProperties` based on the given parameterization and process. This method should have
dispatches added for each process or coupled processes that require additional properties to be defined.
"""
SoilProperties(para::SoilParameterization, proc::Process) = error("not implemented for types $(typeof(para)) and $(typeof(proc))")

"""
    Soil{Tpara<:SoilParameterization,Tprop,Tsp,TP} <: SubSurface{TP}

Generic Soil layer.
"""
Base.@kwdef struct Soil{Tpara<:SoilParameterization,Tprop,Tsp,TP} <: SubSurface{TP}
    proc::TP # subsurface process(es)
    para::Tpara # soil parameterization
    prop::Tprop # soil properties
    sp::Tsp # user-defined specialization
end
"""
    Soil(
        proc::Process;
        para::SoilParameterization=HomogeneousMixture(),
        prop::SoilProperties=SoilProperties(para, proc),
        sp=nothing,
    )

Constructs a `Soil` layer with the given process(es) `proc`, parameterization `para`, and soil properties `prop`.
"""
Soil(
    proc::Process;
    para::SoilParameterization=HomogeneousMixture(),
    prop::SoilProperties=SoilProperties(para, proc),
    sp=nothing,
) = Soil(proc, para, prop, sp)