module Heat

import ConstructionBase

using CryoGrid
using CryoGrid.Hydrology
using CryoGrid.Numerics
using CryoGrid.Utils

using Base: @propagate_inbounds
using FreezeCurves
using FreezeCurves.Solvers

export FreeWater, FreezeCurve

export HeatBalance, ThermalProperties
include("heat_types.jl")

export TemperatureProfile, thermalproperties, freezethaw!, enthalpy, enthalpyinv
include("heat_methods.jl")

export thermalconductivity, thermalconductivity!
include("thermcond.jl")

export heatcapacity, heatcapacity!
include("heatcapacity.jl")

export HeatBC, ConstantTemperature, GeothermalHeatFlux, TemperatureGradient, NFactor
include("heat_bc.jl")

export heatconduction!
include("heat_conduction.jl")

export HeatBalanceImplicit
include("heat_implicit.jl")

include("heat_water.jl")

export StefanProblem, StefanParameters
include("analytic/stefan_analytic.jl")

export heat_conduction_linear_periodic_ub
include("analytic/heat_linear_analytic.jl")

end
