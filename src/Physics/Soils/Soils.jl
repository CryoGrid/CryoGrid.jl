module Soils

using CryoGrid
using CryoGrid.Numerics
using CryoGrid.Utils

using IfElse
using IntervalSets
using ForwardDiff
using FreezeCurves
using FreezeCurves.Solvers
using ModelParameters
using Setfield
using StaticArrays
using Unitful
using UnPack

import CryoGrid
import CryoGrid.InputOutput
import CryoGrid.Heat
import CryoGrid.Hydrology

# from FreezeCurves.jl
export SFCC, PainterKarra, DallAmico, DallAmicoSalt, Westermann, McKenzie
export VanGenuchten, BrooksCorey

# aliases for heat formulations in Heat module
const Temperature = Heat.Temperature
const Enthalpy = Heat.Enthalpy
const EnthalpyImplicit = Heat.EnthalpyImplicit

export Soil, HomogeneousSoil, SoilParameterization
include("types.jl")

export SoilProfile, SoilProperties, soilproperties, porosity, mineral, organic
include("methods.jl")

export HomogeneousMixture, MineralSediment, soilcomponent
include("soil_para.jl")

export RichardsEq
include("soil_water.jl")

export SoilThermalProperties
include("soil_heat.jl")

include("soil_water_heat_coupled.jl")

"""
    Soil(
        proc::Process;
        para::HomogeneousMixture=HomogeneousMixture(),
        solver=default_sfccsolver(proc),
        sp=nothing,
        prop_kwargs...,
    )

Constructs a `HomogeneousSoil` layer with the given process(es) `proc` and parameterization `para`. Additional
keyword arguments are passed through to `SoilProperties`.
"""
Soil(
    proc::Process;
    para::HomogeneousMixture=HomogeneousMixture(),
    solver=default_sfccsolver(proc),
    sp=nothing,
    prop_kwargs...,
) = HomogeneousSoil(para, soilproperties(para, proc; prop_kwargs...), proc, solver, sp)

end
