module Hydrology

import CryoGrid
import ConstructionBase

using CryoGrid
using CryoGrid.Numerics
using CryoGrid.Utils

export WaterBalanceProperties, HydraulicProperties

export WaterBalance, WaterFlow, BucketScheme, NoFlow, Evapotranspiration
include("water_types.jl")

export SaturationProfile,  hydraulicproperties, hydraulicconductivity!
include("water_methods.jl")

include("water_balance.jl")

export DampedET, EvapOnly
include("water_ET.jl")

export ConstantInfiltration, ImpermeableBoundary, WaterBC
include("water_bc.jl")

end
