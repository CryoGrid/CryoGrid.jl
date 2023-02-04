module Hydrology

import CryoGrid
import ConstructionBase

using CryoGrid
using CryoGrid.Physics
using CryoGrid.Numerics
using CryoGrid.Utils

using IfElse
using ModelParameters
using Unitful

export WaterBalanceProperties, HydraulicProperties

export WaterBalance, BucketScheme, NoFlow, Evapotranspiration
include("types.jl")

export SaturationProfile,  hydraulicproperties, hydraulicconductivity!
include("methods.jl")

include("water_balance.jl")

export DampedET, EvapOnly
include("water_ET.jl")

export ConstantInfiltration, Rainfall, WaterBC
include("water_bc.jl")

end
