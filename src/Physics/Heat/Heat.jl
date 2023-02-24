module Heat

import CryoGrid
import CryoGrid.Hydrology

import ConstructionBase

using CryoGrid
using CryoGrid.Numerics
using CryoGrid.Utils

using Base: @propagate_inbounds
using FreezeCurves

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
