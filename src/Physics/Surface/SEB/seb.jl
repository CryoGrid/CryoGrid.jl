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
    solver::TSolver = SimpleNewtonRaphson()
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
    soil = SurfaceProperties(
        α = param(0.2, desc="Initial snow-free surface albedo"),
        ϵ = param(0.97, desc="Initial snow-free surface emissivity"),
        z₀ = param(1e-3, units=u"m", desc="Initial snow-free surface roughness length"),
        rₛ = param(50.0, units=u"s/m", desc="Initial snow-free surface resistance against evapotranspiration and sublimation"),
    ),
    
    snow = SurfaceProperties(
        α = param(0.8, desc="Initial snow-covered surface albedo"),
        ϵ = param(0.99, desc="Initial snow-covered surface emissivity"),
        z₀ = param(5e-4, desc="Initial snow-covered surface roughness length"),
        rₛ = param(50.0, units=u"s/m", desc="Initial snow-covered surface resistance against evapotranspiration and sublimation"),
    ),

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

    z = 2.0u"m"                     # height in meters of air and wind inputs
)

"""
    SurfaceEnergyBalance{TSolution,TStabFun,TPara,TInputs} <: BoundaryProcess{HeatBalance}

Surface energy balance upper boundary condition.
"""
struct SurfaceEnergyBalance{TSolution<:SolutionScheme,TStabFun<:StabilityFunctions,TPara,TInputs} <: BoundaryProcess{HeatBalance}
    inputs::TInputs
    para::TPara
    # type-dependent parameters
    solscheme::TSolution
    stabfun::TStabFun
end
# convenience constructor for SEB
function SurfaceEnergyBalance(;
    Tair=Input(:Tair),
    Lin=Input(:Lin),
    Sin=Input(:Sin),
    pr=Input(:pr),
    qh=Input(:qh),
    wind=Input(:wind),
    para::SEBParams = SEBParams(),
    solscheme::SolutionScheme = Numerical(),
    stabfun::StabilityFunctions = HøgstrømSHEBA(),
)
    inputs = (; Tair, Lin, Sin, pr, qh, wind)
    SurfaceEnergyBalance(inputs, para, solscheme, stabfun)
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
