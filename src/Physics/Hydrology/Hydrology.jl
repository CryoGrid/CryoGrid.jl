module Hydrology

using CryoGrid
using CryoGrid.Numerics
using CryoGrid.Utils

import ConstructionBase
import Interpolations as Interp

export WaterBalanceProperties, HydraulicProperties

export WaterBalance, WaterFlow, BucketScheme, NoFlow, Evapotranspiration
include("water_types.jl")

export ConstantInfiltration, ImpermeableBoundary, WaterBC
include("water_bc.jl")

export hydraulicproperties, hydraulicconductivity!, watercontent!, watercontent, maxwater, minwater, waterdensity
include("water_methods.jl")

export WaterTableInitializer
include("water_init.jl")

include("water_balance.jl")

export DampedET, EvapOnly
include("water_ET.jl")

end
