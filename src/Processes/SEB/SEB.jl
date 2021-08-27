module SEB

using ..HeatConduction: Heat
using ..Processes
using CryoGrid.Forcings
using CryoGrid.Layers: Soil
using CryoGrid.Numerics
using CryoGrid.Utils

using Parameters

import CryoGrid: BoundaryProcess, BoundaryStyle, Neumann, Top
import CryoGrid: initialcondition!, variables

@with_kw struct SEBParams <: Params
    # surface properties --> should be associated with the Stratigraphy and maybe made state variables
    α::Float"1" = 0.2xu"1"                          # surface albedo [-]
    ϵ::Float"1" = 0.97xu"1"                         # surface emissivity [-]
    z₀::Float"m" = 1e-3xu"m"                        # surface roughness length [m]
    rₛ::Float"1/m" = 50.0xu"s/m"                    # surface resistance against evapotranspiration and sublimation [s/m]

    # "natural" constant
    σ::Float"J/(s*m^2*K^4)" = 5.6704e-8xu"J/(s*m^2*K^4)"   # Stefan-Boltzmann constant
    κ::Float"1" = 0.4xu"1"                          # von Kármán constant [-]
    γ::Float"1" = 0.622xu"1"                        # Psychrometric constant [-]
    Rₐ::Float"J/(kg*K)" = 287.058xu"J/(kg*K)"       # specific gas constant of air [J/(kg*K)]
    g::Float"m/s^2" = 9.81xu"m/s^2"                 # gravitational acceleration [m/s^2]

    # material properties (assumed to be constant)
    ρₐ::Float"kg/m^3" = 1.293xu"kg/m^3"             # density of air at standard pressure and 0°C [kg/m^3]
    cₐ::Float"J/(m^3*K)" = 1005.7xu"J/(kg*K)" * ρₐ  # volumetric heat capacity of dry air at standard pressure and 0°C [J/(m^3*K)]
end

struct SurfaceEnergyBalance{F} <: BoundaryProcess
    forcing::F
    sebparams::SEBParams
    function SurfaceEnergyBalance(
        Tair::Forcing,
        p::Forcing,
        q::Forcing,
        wind::Forcing,
        Lin::Forcing,
        Sin::Forcing,
        z::Float"m",
        params::SEBParams=SEBParams()
    )
        forcing = (Tair = Tair, p = p, q = q, wind = wind, Lin = Lin, Sin = Sin, z = z);
        new{typeof(forcing)}(forcing, params)
    end
end

BoundaryStyle(::Type{<:SurfaceEnergyBalance}) = Neumann()

include("seb_heat.jl")

end