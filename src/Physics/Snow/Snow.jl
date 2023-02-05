module Snow

using CryoGrid
using CryoGrid: ContinuousEvent, Increasing, Decreasing # for events/callbacks
using CryoGrid.InputOutput: Forcing
using CryoGrid.Physics
using CryoGrid.Physics.Heat
using CryoGrid.Physics.Hydrology
using CryoGrid.Numerics
using CryoGrid.Utils

import CryoGrid
import CryoGrid.InputOutput
import CryoGrid.Physics
import CryoGrid.Physics.Heat

using IfElse
using ModelParameters
using Unitful
using UnPack

SnowThermalProperties = Heat.ThermalProperties

export Snowpack, SnowProperties, SnowMassBalance
include("types.jl")

include("methods.jl")

export Snowfall
include("snow_bc.jl")

include("snow_bulk.jl")

# Default method implementations for Snowpack type

CryoGrid.basevariables(::Snowpack, ::SnowMassBalance) = (
    Diagnostic(:dsn, Scalar, u"m", domain=0..Inf),
    Diagnostic(:T_ub, Scalar, u"°C"),
)
# for prescribed snow depth/density, the mass balance is given so we do not need to do anything here
CryoGrid.prognosticstep!(::Snowpack, ::SnowMassBalance{<:PrescribedSnow}, ssnow) = nothing
# thermal properties snowpack
Heat.thermalproperties(snow::Snowpack) = snow.prop.heat
# volumetric fractions for snowpack
@inline function Physics.volumetricfractions(::Snowpack, state, i)
    @inbounds let θwi = state.θwi[i],
        θw = state.θw[i],
        θa = 1.0 - θwi,
        θi = θwi - θw;
        return (θw, θi, θa)
    end
end

end