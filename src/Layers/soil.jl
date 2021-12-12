"""
Represents the texture classification of the soil. Sand, Silt, and Clay are provided by default.
"""
abstract type SoilTexture end
struct Sand <: SoilTexture end
struct Silt <: SoilTexture end
struct Clay <: SoilTexture end
"""
    SoilParameterization

Abstract base type for parameterizations of soil properties.
"""
abstract type SoilParameterization end
"""
    SoilCharacteristicFractions{P1,P2,P3,P4} <: SoilParameterization

Represents the composition of the soil in terms of fractions: excess ice, natural porosity, saturation, and organic/(mineral + organic).
"""
struct SoilCharacteristicFractions{P1,P2,P3,P4} <: SoilParameterization
    xic::P1 # excess ice fraction
    por::P2 # natural porosity
    sat::P3 # saturation
    org::P4 # organic fraction of solid; mineral fraction is 1-org
    SoilCharacteristicFractions(xic::P1, por::P2, sat::P3, org::P4) where {P1,P2,P3,P4} = new{P1,P2,P3,P4}(xic,por,sat,org)
end
function soilparameters(;xic=0.0, por=0.5, sat=1.0, org=0.5)
    params = Tuple(Param(p, bounds=(0.0,1.0)) for p in [xic,por,sat,org])
    SoilCharacteristicFractions(params...)
end
SoilProfile(pairs::Pair{<:DistQuantity,<:SoilParameterization}...) = Profile(pairs...)
# Helper functions for obtaining soil compositions from characteristic fractions.
soilcomp(::Val{var}, fracs::SoilCharacteristicFractions) where var = soilcomp(Val{var}(), fracs.xic, fracs.por, fracs.sat, fracs.org)
soilcomp(::Val{:θp}, χ, ϕ, θ, ω) = (1-χ)*ϕ
soilcomp(::Val{:θw}, χ, ϕ, θ, ω) = χ + (1-χ)*ϕ*θ
soilcomp(::Val{:θm}, χ, ϕ, θ, ω) = (1-χ)*(1-ϕ)*(1-ω)
soilcomp(::Val{:θo}, χ, ϕ, θ, ω) = (1-χ)*(1-ϕ)*ω
"""
Thermal conductivity constants.
"""
@with_kw struct SoilTCParams @deftype Float"W/(m*K)"
    kw = 0.57xu"W/(m*K)" #water [Hillel(1982)]
    ko = 0.25xu"W/(m*K)" #organic [Hillel(1982)]
    km = 3.8xu"W/(m*K)" #mineral [Hillel(1982)]
    ka = 0.025xu"W/(m*K)" #air [Hillel(1982)]
    ki = 2.2xu"W/(m*K)" #ice [Hillel(1982)]
end
"""
Heat capacity constants.
"""
@with_kw struct SoilHCParams @deftype Float"J/(K*m^3)"
    cw = 4.2e6xu"J/(K*m^3)" #[J/m^3K] heat capacity water
    co = 2.5e6xu"J/(K*m^3)" #[J/m^3K]  heat capacity organic
    cm = 2e6xu"J/(K*m^3)" #[J/m^3K]  heat capacity mineral
    ca = 0.00125e6xu"J/(K*m^3)" #[J/m^3K]  heat capacity pore space
    ci = 1.9e6xu"J/(K*m^3)" #[J/m^3K]  heat capacity ice
end
"""
Basic Soil layer.
"""
@with_kw struct Soil{T,P<:SoilParameterization,S} <: SubSurface
    texture::T = Sand()
    para::P = soilparameters()
    tc::SoilTCParams = SoilTCParams()
    hc::SoilHCParams = SoilHCParams()
    sp::S = nothing # user-defined specialization
end
# Methods
porosity(soil::Soil{T,<:SoilCharacteristicFractions}) where T = soilcomp(Val{:θp}(), soil.para)
variables(::Soil) = (
    Diagnostic(:θw, Float"1", OnGrid(Cells)), # total water content (xice + saturated pore space)
    Diagnostic(:θm, Float"1", OnGrid(Cells)), # mineral content
    Diagnostic(:θo, Float"1", OnGrid(Cells)), # organic content
)
function initialcondition!(soil::Soil{T,<:SoilCharacteristicFractions}, state) where T
    state.θw .= soilcomp(Val{:θw}(), soil.para)
    state.θm .= soilcomp(Val{:θm}(), soil.para)
    state.θo .= soilcomp(Val{:θo}(), soil.para)
end
