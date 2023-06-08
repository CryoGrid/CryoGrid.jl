module Heat

using CryoGrid
using CryoGrid.Hydrology
using CryoGrid.Numerics
using CryoGrid.Utils

using Base: @propagate_inbounds
using Reexport: @reexport

import ConstructionBase
import Interpolations as Interp

@reexport using FreezeCurves

export FreeWater, FreezeCurve

export HeatBalance, ThermalProperties
include("heat_types.jl")

export TemperatureProfile, thermalproperties, freezethaw!, enthalpy, enthalpyinv
include("heat_methods.jl")

export thermalconductivity, thermalconductivity!
include("thermcond.jl")

export heatcapacity, heatcapacity!
include("heatcapacity.jl")

export HeatBC, ConstantTemperature, GeothermalHeatFlux, TemperatureGradient, GroundHeatFlux, NFactor
include("heat_bc.jl")

export LinearTwoPhaseInitialTempProfile
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
