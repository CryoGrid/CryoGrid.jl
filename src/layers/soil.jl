"""
Represents the general type of the soil. Sand, Silt, and Clay are provided by default.
"""
abstract type SoilType end
struct Sand <: SoilType end
struct Silt <: SoilType end
struct Clay <: SoilType end
"""
Represents the composition type, homogenous or heterogeneous.
"""
abstract type SoilComposition end
struct Homogeneous <: SoilComposition end
struct Heterogeneous <: SoilComposition end
"""
Thermal conductivity constants.
"""
@with_kw struct SoilTCParams <: Params @deftype Float"W/(m*K)"
    kw = 0.57xu"W/(m*K)" #water [Hillel(1982)]
    ko = 0.25xu"W/(m*K)" #organic [Hillel(1982)]
    km = 3.8xu"W/(m*K)" #mineral [Hillel(1982)]
    ka = 0.025xu"W/(m*K)" #air [Hillel(1982)]
    ki = 2.2xu"W/(m*K)" #ice [Hillel(1982)]
end
"""
Heat capacity constants.
"""
@with_kw struct SoilHCParams <: Params @deftype Float"J/(K*m^3)"
    cw = 4.2*10^6xu"J/(K*m^3)" #[J/m^3K] heat capacity water
    co = 2.5*10^6xu"J/(K*m^3)" #[J/m^3K]  heat capacity organic
    cm = 2*10^6xu"J/(K*m^3)" #[J/m^3K]  heat capacity mineral
    ca = 0.00125*10^6xu"J/(K*m^3)" #[J/m^3K]  heat capacity pore space
    ci = 1.9*10^6xu"J/(K*m^3)" #[J/m^3K]  heat capacity ice
end
"""
Parameter type for Soil layers, includes thermal conductivity and heat capacity
constants as well as type/composition.
"""
@with_kw struct SoilParams{TType<:SoilType} <: Params
    tc::SoilTCParams = SoilTCParams()
    hc::SoilHCParams = SoilHCParams()
    type::TType = Sand()
end

function SoilProfile(pairs::Pair{<:DistQuantity,NTuple{5,Float64}}...)
    # order: water+ice (total), liquidWater, mineral, organic, porosity
    @assert begin
        all([(p[2][1] + p[2][3] + p[2][4] ≈ 1.0) || (p[2][5] + p[2][3] + p[2][4] ≈ 1.0) for p in pairs])
    end "either (waterIce + mineral + organic == 1.0) or (porosity + mineral + organic == 1.0) must hold"
    Profile(pairs...;names=(:θw,:θl,:θm,:θo,:θp))
end

"""
Basic Soil layer.
"""
struct Soil{TType<:SoilType,TComp<:SoilComposition} <: SubSurface
    profile::DimArray
    params::SoilParams{TType}
    function Soil(profile::DimArray; kwargs...)
        params = SoilParams(kwargs...)
        shape = size(profile)
        TComp = length(shape) == 1 || shape[1] == 1 ? Homogeneous : Heterogeneous
        if TComp == Homogeneous && length(shape) > 1
            # select depth axis
            profile = profile[Z(1)]
        end
        new{typeof(params.type),TComp}(profile, params)
    end
end

export Soil, SoilProfile, SoilParams
export SoilType, Sand, Silt, Clay
export SoilComposition, Homogeneous, Heterogeneous

Base.show(io::IO, soil::Soil{T}) where T = print(io, "Soil($(soil.params))")

variables(::Soil{T,<:Heterogeneous}) where T = (
    Diagnostic(:θw, Float64, OnGrid(Cells)),
    Diagnostic(:θl, Float64, OnGrid(Cells)),
    Diagnostic(:θm, Float64, OnGrid(Cells)),
    Diagnostic(:θo, Float64, OnGrid(Cells)),
    Diagnostic(:θp, Float64, OnGrid(Cells))
)

variables(soil::Soil{T,<:Homogeneous}) where T = (
    Diagnostic(:θw, Float64, OnGrid(Cells)),
    Diagnostic(:θl, Float64, OnGrid(Cells)),
    Diagnostic(:θm, Float64, OnGrid(Cells)),
    Diagnostic(:θo, Float64, OnGrid(Cells)),
    Diagnostic(:θp, Float64, OnGrid(Cells)),
    Parameter(:θwᵢ, soil.profile[Y(:θw)]),
    Parameter(:θmᵢ, soil.profile[Y(:θm)]),
    Parameter(:θoᵢ, soil.profile[Y(:θo)]),
    Parameter(:θpᵢ, soil.profile[Y(:θp)])
)

function initialcondition!(soil::Soil{T,<:Heterogeneous}, state) where T
    let profile = soil.profile,
        (depths,names) = dims(profile),
        z = ustrip.(depths);
        for p in names
            # in case state is unit-free, reinterpret to match eltype of profile
            pstate = DimArray(similar(state[p], Union{Missing,eltype(profile)}), (Z(state.grids[p]),))
            pstate .= missing
            # assign points where profile is defined
            knots = @view pstate[Z(Near(z))]
            knots .= profile[Z(:),Y(p)]
            # forward fill between points
            state[p] .= Impute.locf(pstate)
        end
    end
end

function initialcondition!(soil::Soil{T,<:Homogeneous}, state) where T
    let profile = soil.profile,
        (names,) = dims(profile);
        for p in names
            if p == :θl
                state[p] .= profile[Y(:θl)]
            else
                state[p] .= state[Symbol(p,:ᵢ)]
            end
        end
    end
end

function thermalconductivity(params::SoilParams, totalWater, liquidWater, mineral, organic)
    @unpack kw,ko,km,ka,ki = params.tc
    let air = 1.0 - totalWater - mineral - organic,
        ice = totalWater - liquidWater,
        water = liquidWater;
        (water*kw^0.5 + ice*ki^0.5 + mineral*km^0.5 + organic*ko^0.5 + air*ka^0.5)^2
    end
end

function heatcapacity(params::SoilParams, totalWater, liquidWater, mineral, organic)
    @unpack cw,co,cm,ca,ci = params.hc
    let air = 1.0 - totalWater - mineral - organic,
        ice = totalWater - liquidWater,
        water = liquidWater;
        water*cw + ice*ci + mineral*cm + organic*co + air*ca
    end
end

export thermalconductivity, heatcapacity
