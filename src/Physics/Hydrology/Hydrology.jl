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

"""
    WaterBalanceProperties

Numerical constants shared across water balance implementations.
"""
Base.@kwdef struct WaterBalanceProperties{Tρw,TLsg,Trb,Trc}
    ρw::Tρw = Physics.Constants.ρw
    Lsg::TLsg = Physics.Constants.Lsg
    r_β::Trb = 1e3 # reduction factor scale parameter
    r_c::Trc = 0.96325 # reduction factor shift parameter
end
CryoGrid.parameterize(prop::WaterBalanceProperties) = prop
"""
    HydraulicProperties

Material hydraulic properties.
"""
Utils.@properties HydraulicProperties(
    kw_sat = 1e-5u"m/s",
)
hydraulicproperties(sub::SubSurface) = error("hydraulic properties not defined for layer $(typeof(sub))")

function CryoGrid.parameterize(prop::HydraulicProperties)
    return HydraulicProperties(
        map(values(prop)) do val
            # this currently assumes that all properties have a strictly positive domain!
            CryoGrid.parameterize(val, domain=StrictlyPositive)
        end
    )
end

export BucketScheme, Evapotranspiration, WaterBalance
include("water_balance.jl")
export DampedET, EvapOnly
include("water_ET.jl")
export ConstantInfiltration, Rainfall, WaterBC
include("water_bc.jl")

end
