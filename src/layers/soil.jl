abstract type SoilType end
struct Sand end
struct Silt end
struct Clay end

@with_kw struct SoilTCParams <: Params @deftype UFloat"W/(m*K)"
    ka = 0.025U"W/(m*K)" #air [Hillel(1982)]
    kw = 0.57U"W/(m*K)" #water [Hillel(1982)]
    ko = 0.25U"W/(m*K)" #organic [Hillel(1982)]
    km = 3.8U"W/(m*K)" #mineral [Hillel(1982)]
    ki = 2.2U"W/(m*K)" #ice [Hillel(1982)]
end
@with_kw struct SoilHCParams <: Params @deftype UFloat"J/(K*m^3)"
    cw = 4.2*10^6U"J/(K*m^3)" #[J/m^3K] heat capacity water
    co = 2.5*10^6U"J/(K*m^3)" #[J/m^3K]  heat capacity organic
    cm = 2*10^6U"J/(K*m^3)", #[J/m^3K]  heat capacity mineral
    ca = 0.00125*10^6U"J/(K*m^3)" #[J/m^3K]  heat capacity pore space
    ci = 1.9*10^6U"J/(K*m^3)" #[J/m^3K]  heat capacity ice
end

const SoilProfile{D,Q} = Profile{D,5,Q} where {D,Q<:DistQuantity}
function SoilProfile(pairs::Pair{Q,NTuple{5,Float64}}...) where {Q<:DistQuantity}
    # order: water+ice (total), liquidWater, mineral, organic
    @assert all([p.second[1] + p.second[3] + p.second[4] <= 1.0 for p in pairs]) "composition must be <= 1.0"
    @assert all([p.second[end] >= 0.0 for p in pairs]) "porosity must be >= 0.0"
    Profile(pairs...;names=(:θ_w,:θ_l,:θ_m,:θ_o,:por))
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
    Var(:θ_w, Float64, OnGrid(Cells)),
    Var(:θ_l, Float64, OnGrid(Cells)),
    Var(:θ_m, Float64, OnGrid(Cells)),
    Var(:θ_o, Float64, OnGrid(Cells)),
    Var(:por, Float64, OnGrid(Cells))
)

function thermalConductivity(params::SoilTCParams, totalWater, mineral, organic)
    @unpack kw, km, ko, ka = params
    air = 1.0 - totalWater - mineral - organic;
    (water*kw^0.5 + ice*ki^0.5 + mineral*km^0.5 + organic*ko^0.5 + air*ka^0.5)^2
end

function heatCapacity(params::SoilHCParams, totalWater, mineral, organic)
    @unpack cw, ci, cm, co, ca = params
    air = 1.0 - waterIce - mineral - organic
    water*cw + ice*ci + mineral*cm + organic*co + air*ca
end

export Soil, SoilProfile, SoilTCParams, SoilHCParams
export SoilType, Sand, Silt, Clay
export thermalConductivity, heatCapacity, diagnostic_vars
