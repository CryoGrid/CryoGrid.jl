module Hydrology

using CryoGrid
using CryoGrid.Numerics
using CryoGrid.Utils

import ConstructionBase
import Interpolations as Interp

export WaterBalanceProperties, HydraulicProperties

export WaterBalance, WaterFlow, BucketScheme, NoFlow, Evapotranspiration
include("water_types.jl")

export SaturationProfile,  hydraulicproperties, hydraulicconductivity!
include("water_methods.jl")

export ConstantInfiltration, ImpermeableBoundary, WaterBC
include("water_bc.jl")

export WaterTableInitializer
include("water_init.jl")

include("water_balance.jl")

export DampedET, EvapOnly
include("water_ET.jl")

end
