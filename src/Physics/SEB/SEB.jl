module SEB

using CryoGrid.InputOutput: Forcing
using CryoGrid.Physics
using CryoGrid.Physics.Boundaries
using CryoGrid.Physics.HeatConduction
using CryoGrid.Physics.Hydrology
using CryoGrid.Physics.Soils
using CryoGrid.Numerics
using CryoGrid.Utils

import CryoGrid: BoundaryProcess, BoundaryStyle, Neumann, Top
import CryoGrid: initialcondition!, variables, boundaryvalue

using Unitful

export SurfaceEnergyBalance, SEBParams
export Businger, HøgstrømSHEBA, Iterative, Analytical, Numerical

abstract type SolutionScheme end
"""
Byun 1990
"""
struct Analytical <: SolutionScheme end

"""
equation by Westermann2016, but use Newton solver
- not implemented yet
"""
struct Numerical <: SolutionScheme end

"""
Westermann2016, use info from last time step
"""
struct Iterative <: SolutionScheme end

abstract type StabilityFunctions end

"""
Businger 1971
"""
struct Businger <: StabilityFunctions end

"""
Høgstrøm 1988 (unstable conditions)
SHEBA, Uttal et al., 2002, Grachev et al. 2007 (stable conditions)
"""
struct HøgstrømSHEBA <: StabilityFunctions end

Base.@kwdef struct SEBParams{TSolution,TStabFun}
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

    Pr₀::Float64 = 0.74                             # turbulent Prandtl number
    βₘ::Float64 = 4.7
    βₕ::Float64 = βₘ/Pr₀
    γₘ::Float64 = 15.0
    γₕ::Float64 = 9.0

    # type-dependent parameters
    solscheme::TSolution = Iterative()
    stabfun::TStabFun = HøgstrømSHEBA()
end

"""
    SurfaceEnergyBalance{TSolution,TStabFun,F} <: BoundaryProcess{Heat}

Surface energy balance upper boundary condition.
"""
struct SurfaceEnergyBalance{TSolution,TStabFun,F} <: BoundaryProcess{Heat}
    forcings::F
    sebparams::SEBParams{TSolution,TStabFun}
    SurfaceEnergyBalance(forcings::NamedTuple, sebparams::SEBParams) = new{typeof(sebparams.solscheme),typeof(sebparams.stabfun),typeof(forcings)}(forcings, sebparams)
    function SurfaceEnergyBalance(
        Tair::Forcing{u"°C"},
        p::Forcing{u"kPa"},
        q::Forcing{u"kg/kg"},
        wind::Forcing{u"m/s"},
        Lin::Forcing{u"W/m^2"},
        Sin::Forcing{u"W/m^2"},
        z;
        kwargs...
    )
        forcings = (; Tair, p, q, wind, Lin, Sin, z);
        sebparams = SEBParams(;kwargs...);
        new{typeof(sebparams.solscheme),typeof(sebparams.stabfun),typeof(forcings)}(forcings, sebparams)
    end
end

include("seb_heat.jl")
include("seb_water.jl")

end
