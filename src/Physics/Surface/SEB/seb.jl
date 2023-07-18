abstract type SolutionScheme end
"""
    Analytical

Analytical formulation from Byun 1990.
"""
struct Analytical <: SolutionScheme end

"""
    Iterative

Westermann 2016, use state from last time step.
"""
struct Iterative <: SolutionScheme end

"""
    Numerical

Equations from Westermann 2016, but use nonlinear solver.
"""
Base.@kwdef struct Numerical{TSolver} <: SolutionScheme
    solver::TSolver = NewtonRaphson()
end
Flatten.flattenable(::Type{<:Numerical}, ::Type) = false

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

Utils.@properties SurfaceProperties(
    # surface properties --> should be associated with the Stratigraphy and maybe made state variables or parameters
    α = 0.2,             # initial surface albedo [-]
    ϵ = 0.97,            # initial surface emissivity [-]
    z₀ = 1e-3u"m",       # initial surface roughness length [m]
    rₛ = 50.0u"s/m",     # initial surface resistance against evapotranspiration and sublimation [s/m]
)

Utils.@properties SEBParams(
    # surface propertiess; includes soil and snow by default;
    # this can be easily extended by the user to include other layers by adding more
    # fields to `SEBParams` and then implementing `surfaceproperties` for those layers.
    soil = SurfaceProperties(),
    snow = SurfaceProperties(α=0.8, ϵ=0.99, z₀=5e-4),

    # "natural" constant
    σ = 5.6704e-8u"J/(s*m^2*K^4)",   # Stefan-Boltzmann constant
    κ = 0.4u"1",                     # von Kármán constant [-]
    γ = 0.622u"1",                   # Psychrometric constant [-]
    Rₐ = 287.058u"J/(kg*K)",         # specific gas constant of air [J/(kg*K)]
    g = 9.81u"m/s^2",                # gravitational acceleration [m/s^2]

    # material properties (assumed to be constant)
    ρₐ = 1.293u"kg/m^3",             # density of air at standard pressure and 0°C [kg/m^3]
    cₐ = 1005.7u"J/(kg*K)" * ρₐ,     # volumetric heat capacity of dry air at standard pressure and 0°C [J/(m^3*K)]

    Pr₀ = 0.74,                      # turbulent Prandtl number
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
    # Default constructor with exact fields for reconstruction
    SurfaceEnergyBalance(forcings::NamedTuple, para::SEBParams, solscheme::SolutionScheme, stabfun::StabilityFunctions) =
        new{typeof(solscheme),typeof(stabfun),typeof(para),typeof(forcings)}(forcings, para, solscheme, stabfun)
    # User facing constructors
    function SurfaceEnergyBalance(
        Tair::TemperatureForcing, # air temperature
        pr::PressureForcing, # air pressure
        qh::HumidityForcing, # specific humidity
        wind::VelocityForcing, # non-directional wind speed
        Lin::EnergyFluxForcing, # long-wave incoming radiation
        Sin::EnergyFluxForcing, # short-wave incoming radiation
        z; # height [m] of air temperature and wind forcing
        para::SEBParams = SEBParams(),
        solscheme::SolutionScheme = Numerical(),
        stabfun::StabilityFunctions = HøgstrømSHEBA(),
    )
        forcings = (; Tair, pr, qh, wind, Lin, Sin, z);
        SurfaceEnergyBalance(forcings, para, solscheme, stabfun)
    end
end

"""
    surfaceproperties(::SurfaceEnergyBalance, ::SubSurface)

Retrieves the `SurfaceProperties` for the given `SubSurface` layer.
"""
surfaceproperties(seb::SurfaceEnergyBalance, sub::SubSurface) = error("surfaceproperties not implemented for layer of type $(typeof(sub))")
surfaceproperties(seb::SurfaceEnergyBalance, ::Soil) = seb.para.soil
surfaceproperties(seb::SurfaceEnergyBalance, ::Snowpack) = seb.para.snow

include("seb_state.jl")
include("seb_solve.jl")
include("seb_heat.jl")
