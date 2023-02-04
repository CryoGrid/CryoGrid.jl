module Heat

using CryoGrid
using CryoGrid.InputOutput: Forcing
using CryoGrid.Physics
using CryoGrid.Physics.Hydrology
using CryoGrid.Numerics
using CryoGrid.Utils

using Base: @propagate_inbounds
using IfElse
using FreezeCurves: FreezeCurves, FreezeCurve, FreeWater
using ModelParameters
using Unitful
using UnPack

import CryoGrid
import CryoGrid.Physics

import ConstructionBase

export FreeWater, FreezeCurve

export HeatBalance, ThermalProperties
include("types.jl")

export TemperatureProfile, thermalproperties, freezethaw!, enthalpy, enthalpyinv
include("methods.jl")

export thermalconductivity, thermalconductivity!
include("thermcond.jl")

export heatcapacity, heatcapacity!
include("heatcap.jl")

export HeatBC, ConstantTemperature, GeothermalHeatFlux, TemperatureGradient, NFactor
include("heat_bc.jl")

export heatconduction!
include("heat_conduction.jl")

include("heat_water.jl")

export HeatBalanceImplicit
include("heat_implicit.jl")

export StefanProblem, StefanParameters
include("analytic/stefan_analytic.jl")

end
