"""
    SoilParameterization

Base type for parameterizations of soil consituents.
"""
abstract type SoilParameterization end

"""
    Soil{Tpara<:SoilParameterization,Tprop,Tsp,TP} <: SubSurface{TP}

Generic Soil layer.
"""
struct Soil{Tpara<:SoilParameterization,Tprop,Tsp,TP} <: SubSurface{TP}
    proc::TP # subsurface process(es)
    para::Tpara # soil parameterization
    prop::Tprop # soil properties
    sp::Tsp # user-defined specialization
end
"""
    Soil(
        proc::Process;
        para::SoilParameterization=HomogeneousMixture(),
        sp=nothing,
        prop_kwargs...,
    )

Constructs a `Soil` layer with the given process(es) `proc` and parameterization `para`. Additional
keyword arguments are passed through to `SoilProperties`.
"""
Soil(
    proc::Process;
    para::SoilParameterization=HomogeneousMixture(),
    sp=nothing,
    prop_kwargs...,
) = Soil(proc, para, soilproperties(para, proc; prop_kwargs...), sp)
