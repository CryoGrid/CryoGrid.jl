module Hydrology

import CryoGrid
import ConstructionBase

using CryoGrid
using CryoGrid.Physics
using CryoGrid.Numerics
using CryoGrid.Numerics: flux!, divergence!, ∇
using CryoGrid.Utils

using IfElse
using ModelParameters
using Unitful

export HydraulicProperties

Base.@kwdef struct HydraulicProperties{Tconsts,Tkwsat,Trb,Trc}
    consts::Tconsts = Physics.Constants()
    kw_sat::Tkwsat = Param(1e-5, domain=0..Inf, units=u"m/s")
    r_β::Trb = 1e3 # reduction factor scale parameter
    r_c::Trc = 0.96325 # reduction factor shift parameter
end

export BucketScheme, WaterBalance
include("water_balance.jl")
export ConstantInfiltration, Rainfall, WaterBC
include("water_bc.jl")

end
