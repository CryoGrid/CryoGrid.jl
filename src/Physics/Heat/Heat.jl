module Heat

using CryoGrid
using CryoGrid.Hydrology
using CryoGrid.Numerics
using CryoGrid.Utils

using Base: @propagate_inbounds
using Reexport: @reexport

using IfElse

import ConstructionBase
import Interpolations as Interp

@reexport using FreezeCurves

export FreeWater, FreezeCurve

export HeatBalance
include("heat_types.jl")

export TemperatureProfile, thermalproperties, freezethaw!, freezecurve, enthalpy, enthalpyinv
include("heat_methods.jl")

export ThermalProperties, thermalconductivity, thermalconductivity!, heatcapacity, heatcapacity!
include("thermal_properties.jl")

export HeatBC, ConstantTemperature, GeothermalHeatFlux, TemperatureBC, GroundHeatFlux, NFactor
include("heat_bc.jl")

export PermafrostTemperatureInit, ThermalSteadyStateInit
include("heat_init.jl")

export heatconduction!
include("heat_conduction.jl")

export HeatBalanceImplicit
include("heat_implicit.jl")

export WaterHeatBC
include("heat_water.jl")

export StefanProblem, StefanParameters
include("analytic/stefan_analytic.jl")

export heat_conduction_linear_periodic_ub
include("analytic/heat_linear_analytic.jl")

end
