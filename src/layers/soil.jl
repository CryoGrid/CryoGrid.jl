abstract type SoilType end
struct Sand <: SoilType end
struct Silt <: SoilType end
struct Clay <: SoilType end

@with_kw struct SoilTCParams <: Params @deftype UFloat"W/(m*K)"
    kw = 0.57U"W/(m*K)" #water [Hillel(1982)]
    ko = 0.25U"W/(m*K)" #organic [Hillel(1982)]
    km = 3.8U"W/(m*K)" #mineral [Hillel(1982)]
    ka = 0.025U"W/(m*K)" #air [Hillel(1982)]
    ki = 2.2U"W/(m*K)" #ice [Hillel(1982)]
end
@with_kw struct SoilHCParams <: Params @deftype UFloat"J/(K*m^3)"
    cw = 4.2*10^6U"J/(K*m^3)" #[J/m^3K] heat capacity water
    co = 2.5*10^6U"J/(K*m^3)" #[J/m^3K]  heat capacity organic
    cm = 2*10^6U"J/(K*m^3)" #[J/m^3K]  heat capacity mineral
    ca = 0.00125*10^6U"J/(K*m^3)" #[J/m^3K]  heat capacity pore space
    ci = 1.9*10^6U"J/(K*m^3)" #[J/m^3K]  heat capacity ice
end

const SoilProfile{D,Q,T} = Profile{D,5,Q,T} where {D,Q,T}
function SoilProfile(pairs::Pair{<:DistQuantity,NTuple{5,Float64}}...)
    # order: water+ice (total), liquidWater, mineral, organic
    @assert all([p.second[1] + p.second[3] + p.second[4] <= 1.0 for p in pairs]) "composition must be <= 1.0"
    @assert all([p.second[end] >= 0.0 for p in pairs]) "porosity must be >= 0.0"
    Profile(pairs...;names=(:θw,:θl,:θm,:θo,:por))
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
    Diagnostic(:por, Float64, OnGrid(Cells))
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
