"""
Represents the texture classification of the soil. Sand, Silt, and Clay are provided by default.
"""
abstract type SoilTexture end
struct Sand <: SoilTexture end
struct Silt <: SoilTexture end
struct Clay <: SoilTexture end
"""
Represents the composition of the soil in terms of fractions: excess ice, natural porosity, saturation, and organic/(mineral + organic).
"""
struct SoilComposition{P}
    xic::P # excess ice fraction
    por::P # natural porosity
    sat::P # saturation
    org::P # organic fraction of solid; mineral fraction is 1-org
    SoilComposition(xic::P, por::P, sat::P, org::P) where P = new{P}(xic,por,sat,org)
    function SoilComposition(;xic=0.0, por=0.5, sat=1.0, org=0.5)
        params = Tuple(Param(p, bounds=(0.0,1.0)) for p in [xic,por,sat,org])
        new{typeof(params[1])}(params...)
    end
end
SoilProfile(pairs::Pair{<:DistQuantity,<:SoilComposition}...) = Profile(pairs...)
# Helper functions for obtaining soil component fractions from soil properties.
soilfrac(::Val{:θp}, χ, ϕ, θ, ω) = (1-χ)*ϕ
soilfrac(::Val{:θw}, χ, ϕ, θ, ω) = χ + (1-χ)*ϕ*θ
soilfrac(::Val{:θm}, χ, ϕ, θ, ω) = (1-χ)*(1-ϕ)*(1-ω)
soilfrac(::Val{:θo}, χ, ϕ, θ, ω) = (1-χ)*(1-ϕ)*ω
θx(comp::SoilComposition) = comp.xic
θp(comp::SoilComposition) = soilfrac(Val{:θp}(), comp.xic, comp.por, comp.sat, comp.org)
θw(comp::SoilComposition) = soilfrac(Val{:θw}(), comp.xic, comp.por, comp.sat, comp.org)
θm(comp::SoilComposition) = soilfrac(Val{:θm}(), comp.xic, comp.por, comp.sat, comp.org)
θo(comp::SoilComposition) = soilfrac(Val{:θo}(), comp.xic, comp.por, comp.sat, comp.org)
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
@with_kw struct Soil{T,P,S} <: SubSurface
    texture::T = Sand()
    comp::SoilComposition{P} = SoilComposition()
    tc::SoilTCParams = SoilTCParams()
    hc::SoilHCParams = SoilHCParams()
    sp::S = nothing # user-defined specialization
end
