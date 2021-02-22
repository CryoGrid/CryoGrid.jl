abstract type SoilType end
struct Sand <: SoilType end
struct Silt <: SoilType end
struct Clay <: SoilType end

@with_kw struct SoilTCParams <: Params @deftype Float"W/(m*K)"
    kw = 0.57xu"W/(m*K)" #water [Hillel(1982)]
    ko = 0.25xu"W/(m*K)" #organic [Hillel(1982)]
    km = 3.8xu"W/(m*K)" #mineral [Hillel(1982)]
    ka = 0.025xu"W/(m*K)" #air [Hillel(1982)]
    ki = 2.2xu"W/(m*K)" #ice [Hillel(1982)]
end
@with_kw struct SoilHCParams <: Params @deftype Float"J/(K*m^3)"
    cw = 4.2*10^6xu"J/(K*m^3)" #[J/m^3K] heat capacity water
    co = 2.5*10^6xu"J/(K*m^3)" #[J/m^3K]  heat capacity organic
    cm = 2*10^6xu"J/(K*m^3)" #[J/m^3K]  heat capacity mineral
    ca = 0.00125*10^6xu"J/(K*m^3)" #[J/m^3K]  heat capacity pore space
    ci = 1.9*10^6xu"J/(K*m^3)" #[J/m^3K]  heat capacity ice
end

const SoilProfile{D,Q,T} = Profile{D,5,Q,T} where {D,Q,T}
function SoilProfile(pairs::Pair{<:DistQuantity,NTuple{5,Float64}}...)
    # order: water+ice (total), liquidWater, mineral, organic, porosity
    @assert begin
        all([(p[2][1] + p[2][3] + p[2][4] ≈ 1.0) || (p[2][5] + p[2][3] + p[2][4] ≈ 1.0) for p in pairs])
    end "either (waterIce + mineral + organic == 1.0) or (porosity + mineral + organic == 1.0) must hold"
    Profile(pairs...;names=(:θw,:θl,:θm,:θo,:θp))
end

"""
Basic Soil layer with interchangeable type T to allow for specificity in dispatch.
"""
struct Soil{T,P} <: SubSurface
    profile::P
    tcparams::SoilTCParams
    hcparams::SoilHCParams
    Soil{T}(profile::P, tcparams::SoilTCParams=SoilTCParams(), hcparams=SoilHCParams()) where
        {T<:SoilType,P<:SoilProfile} = new{T,P}(profile,tcparams,hcparams)
end

variables(soil::Soil) = (
    Diagnostic(:θw, Float64, OnGrid(Cells)),
    Diagnostic(:θl, Float64, OnGrid(Cells)),
    Diagnostic(:θm, Float64, OnGrid(Cells)),
    Diagnostic(:θo, Float64, OnGrid(Cells)),
    Diagnostic(:θp, Float64, OnGrid(Cells))
)

function initialcondition!(soil::Soil, state)
    interpolateprofile!(soil.profile, state)
end

function thermalConductivity(params::SoilTCParams, totalWater, liquidWater, mineral, organic)
    @unpack kw, ko, km, ka, ki = params
    let air = 1.0 - totalWater - mineral - organic,
        ice = totalWater - liquidWater,
        water = liquidWater;
        (water*kw^0.5 + ice*ki^0.5 + mineral*km^0.5 + organic*ko^0.5 + air*ka^0.5)^2
    end
end

function heatCapacity(params::SoilHCParams, totalWater, liquidWater, mineral, organic)
    @unpack cw, co, cm, ca, ci = params
    let air = 1.0 - totalWater - mineral - organic,
        ice = totalWater - liquidWater,
        water = liquidWater;
        water*cw + ice*ci + mineral*cm + organic*co + air*ca
    end
end

export Soil, SoilProfile, SoilTCParams, SoilHCParams
export SoilType, Sand, Silt, Clay
export thermalConductivity, heatCapacity, initialcondition!, variables
