module SEB

using CryoGrid
using CryoGrid.Physics
using CryoGrid.Physics
using CryoGrid.Physics.Heat
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

Utils.@properties SEBParams(
    # surface properties --> should be associated with the Stratigraphy and maybe made state variables or parameters
    α = 0.2,             # surface albedo [-]
    ϵ = 0.97,            # surface emissivity [-]
    z₀ = 1e-3u"m",       # surface roughness length [m]
    rₛ = 50.0u"s/m",     # surface resistance against evapotranspiration and sublimation [s/m]

    # "natural" constant
    σ = 5.6704e-8u"J/(s*m^2*K^4)",   # Stefan-Boltzmann constant
    κ = 0.4u"1",                          # von Kármán constant [-]
    γ = 0.622u"1",                        # Psychrometric constant [-]
    Rₐ = 287.058u"J/(kg*K)",       # specific gas constant of air [J/(kg*K)]
    g = 9.81u"m/s^2",                 # gravitational acceleration [m/s^2]

    # material properties (assumed to be constant)
    ρₐ = 1.293u"kg/m^3",             # density of air at standard pressure and 0°C [kg/m^3]
    cₐ = 1005.7u"J/(kg*K)" * ρₐ,     # volumetric heat capacity of dry air at standard pressure and 0°C [J/(m^3*K)]

    Pr₀ = 0.74,                             # turbulent Prandtl number
    βₘ = 4.7,
    βₕ = βₘ/Pr₀,
    γₘ = 15.0,
    γₕ = 9.0,
)

"""
    SurfaceEnergyBalance{TSolution,TStabFun,TPara,F} <: BoundaryProcess{HeatBalance}

Surface energy balance upper boundary condition.
"""
struct SurfaceEnergyBalance{TSolution,TStabFun,TPara,F} <: BoundaryProcess{HeatBalance}
    forcings::F
    para::TPara
    # type-dependent parameters
    solscheme::TSolution
    stabfun::TStabFun
    SurfaceEnergyBalance(forcings::NamedTuple, para::SEBParams, solscheme::SolutionScheme, stabfun::StabilityFunctions) =
        new{typeof(solscheme),typeof(stabfun),typeof(para),typeof(forcings)}(forcings, para, solscheme, stabfun)
    function SurfaceEnergyBalance(
        Tair::Forcing{u"°C"},
        p::Forcing{u"kPa"},
        q::Forcing{u"kg/kg"},
        wind::Forcing{u"m/s"},
        Lin::Forcing{u"W/m^2"},
        Sin::Forcing{u"W/m^2"},
        z;
        para::SEBParams = SEBParams(),
        solscheme::SolutionScheme = Iterative(),
        stabfun::StabilityFunctions = HøgstrømSHEBA(),
    )
        forcings = (; Tair, p, q, wind, Lin, Sin, z);
        SurfaceEnergyBalance(forcings, para, solscheme, stabfun)
    end
end

include("seb_heat.jl")
include("seb_water.jl")

end
