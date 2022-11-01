module Soils

using CryoGrid
using CryoGrid.Numerics
using CryoGrid.Physics
using CryoGrid.Physics.Heat
using CryoGrid.Physics.Hydrology
using CryoGrid.Utils

using IfElse
using IntervalSets
using ForwardDiff
using FreezeCurves
using ModelParameters
using Setfield
using StaticArrays
using Unitful
using UnPack

import CryoGrid
import CryoGrid.InputOutput
import CryoGrid.Physics
import CryoGrid.Physics.Heat
import CryoGrid.Physics.Hydrology

export Soil, SoilParameterization, SoilProperties
export soilcomponent, porosity, mineral, organic

# from FreezeCurves
export SFCC, PainterKarra, DallAmico, DallAmicoSalt, Westermann, McKenzie
export VanGenuchten, BrooksCorey

# aliases for heat formulations in Heat module
const Temperature = Heat.Temperature
const Enthalpy = Heat.Enthalpy
const EnthalpyImplicit = Heat.EnthalpyImplicit

"""
    SoilParameterization

Base type for parameterizations of soil consituents.
"""
abstract type SoilParameterization end

export HomogeneousMixture, SoilProfile
include("soilpara.jl")

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
    Soil{Tpara<:SoilParameterization,Tprop,Tsp,Tproc} <: SubSurface{Tproc}

Generic Soil layer.
"""
Base.@kwdef struct Soil{Tpara<:SoilParameterization,Tprop,Tsp,Tproc} <: SubSurface{Tproc}
    proc::Tproc # subsurface process(es)
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
"""
    mineral(soil::Soil, state, i)

Retrieves the mineral content for the given layer at grid cell `i`, if provided.
Defaults to using the scalar mineral content defined on `soil`.
"""
mineral(soil::Soil, state, i) = Utils.getscalar(mineral(soil, state), i)
"""
    organic(soil::Soil, state, i)

Retrieves the organic content for the given layer at grid cell `i`, if provided.
Defaults to using the scalar organic content defined on `soil`.
"""
organic(soil::Soil, state, i) = Utils.getscalar(organic(soil, state), i)
"""
    porosity(soil::Soil, state, i)

Retrieves the porosity for the given layer at grid cell `i`, if provided.
Defaults to using the scalar porosity defined on `soil`.
"""
porosity(soil::Soil, state, i) = Utils.getscalar(porosity(soil, state), i)

## Homogeneous soils
porosity(soil::Soil{<:HomogeneousMixture}, state=nothing) = soilcomponent(Val{:θp}(), soil.para)
mineral(soil::Soil{<:HomogeneousMixture}, state=nothing) = soilcomponent(Val{:θm}(), soil.para)
organic(soil::Soil{<:HomogeneousMixture}, state=nothing) = soilcomponent(Val{:θo}(), soil.para)
## Heterogeneous soils  
porosity(::Soil, state) = state.θp
mineral(::Soil, state) = state.θm
organic(::Soil, state) = state.θo

CryoGrid.variables(soil::Soil{<:HomogeneousMixture}) = (
    Diagnostic(:θwi, OnGrid(Cells), domain=0..1),
    CryoGrid.variables(soil, processes(soil))...
)
CryoGrid.variables(soil::Soil) = (
    Diagnostic(:θwi, OnGrid(Cells), domain=0..1),
    Diagnostic(:θp, OnGrid(Cells), domain=0..1),
    Diagnostic(:θm, OnGrid(Cells), domain=0..1),
    Diagnostic(:θo, OnGrid(Cells), domain=0..1),
    CryoGrid.variables(soil, processes(soil))...,
)
function CryoGrid.initialcondition!(soil::Soil{<:HomogeneousMixture}, state)
    ϕ = soil.para.por
    θ = soil.para.sat
    ω = soil.para.org
    @. state.θwi = soilcomponent(Val{:θwi}(), ϕ, θ, ω)
    CryoGrid.initialcondition!(soil, processes(soil), state)
end
function CryoGrid.initialcondition!(soil::Soil{<:MineralSediment}, state)
    evaluate(x::Number) = x
    evaluate(f::Function) = f(soil, state)
    θp = evaluate(soil.para.por)
    @. state.θwi = θp # assuming saturated conditions
    @. state.θp = θp
    @. state.θm = 1 - state.θp
    @. state.θo = zero(eltype(state.θo))
    CryoGrid.initialcondition!(soil, processes(soil), state)
end

"""
    SoilProfile(pairs::Pair{<:DistQuantity,<:SoilParameterization}...)

Alias for `Profile(pairs...)` assigning soil parameterizations to specific depths.
"""
SoilProfile(pairs::Pair{<:DistQuantity,<:SoilParameterization}...) = Profile(pairs...)

export RichardsEq
include("soil_water.jl")
include("soil_heat.jl")
include("soil_water_heat_coupled.jl")

end
