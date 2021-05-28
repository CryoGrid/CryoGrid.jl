abstract type SoilType end
struct Sand <: SoilType end
struct Silt <: SoilType end
struct Clay <: SoilType end

abstract type SoilComposition end
struct Homogeneous <: SoilComposition end
struct Heterogeneous <: SoilComposition end

@with_kw struct SoilParams{TComp<:SoilComposition,TType<:SoilType} <: Params
    kw::Float"W/(m*K)" = 0.57xu"W/(m*K)" #water [Hillel(1982)]
    ko::Float"W/(m*K)" = 0.25xu"W/(m*K)" #organic [Hillel(1982)]
    km::Float"W/(m*K)" = 3.8xu"W/(m*K)" #mineral [Hillel(1982)]
    ka::Float"W/(m*K)" = 0.025xu"W/(m*K)" #air [Hillel(1982)]
    ki::Float"W/(m*K)" = 2.2xu"W/(m*K)" #ice [Hillel(1982)]
    cw::Float"J/(K*m^3)" = 4.2*10^6xu"J/(K*m^3)" #[J/m^3K] heat capacity water
    co::Float"J/(K*m^3)" = 2.5*10^6xu"J/(K*m^3)" #[J/m^3K]  heat capacity organic
    cm::Float"J/(K*m^3)" = 2*10^6xu"J/(K*m^3)" #[J/m^3K]  heat capacity mineral
    ca::Float"J/(K*m^3)" = 0.00125*10^6xu"J/(K*m^3)" #[J/m^3K]  heat capacity pore space
    ci::Float"J/(K*m^3)" = 1.9*10^6xu"J/(K*m^3)" #[J/m^3K]  heat capacity ice
    type::TType = Sand()
    comp::TComp = Heterogeneous()
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
struct Soil{TComp<:SoilComposition,TType<:SoilType} <: SubSurface
    profile::DimArray
    params::SoilParams{TComp,TType}
    function Soil(profile::DimArray; kwargs...)
        params = SoilParams(kwargs...)
        @assert typeof(params.type) <: Heterogeneous || size(profile,1) == 1 "Homogeneous soil profile must have exactly one depth/row."
        new{typeof(params.type),typeof(params.comp)}(profile, params)
    end
end

export Soil, SoilProfile, SoilParams
export SoilType, Sand, Silt, Clay
export SoilComposition, Homogeneous, Heterogeneous

Base.show(io::IO, soil::Soil{T}) where T = print(io, "Soil($(soil.params))")

variables(soil::Soil{<:Heterogeneous}) = (
    Diagnostic(:θw, Float64, OnGrid(Cells)),
    Diagnostic(:θl, Float64, OnGrid(Cells)),
    Diagnostic(:θm, Float64, OnGrid(Cells)),
    Diagnostic(:θo, Float64, OnGrid(Cells)),
    Diagnostic(:θp, Float64, OnGrid(Cells))
)

variables(soil::Soil{<:Homogeneous}) = (
    Diagnostic(:θl, Float64, OnGrid(Cells)),
    Parameter(:θw, soil.profile[Z(1),Y(:θw)]),
    Parameter(:θm, soil.profile[Z(1),Y(:θm)]),
    Parameter(:θo, soil.profile[Z(1),Y(:θo)]),
    Parameter(:θp, soil.profile[Z(1),Y(:θp)]),
)

function initialcondition!(soil::Soil{<:Heterogeneous}, state)
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

function thermalconductivity(params::SoilParams, totalWater, liquidWater, mineral, organic)
    @unpack kw,ko,km,ka,ki,_,_,_,_,_ = params
    let air = 1.0 - totalWater - mineral - organic,
        ice = totalWater - liquidWater,
        water = liquidWater;
        (water*kw^0.5 + ice*ki^0.5 + mineral*km^0.5 + organic*ko^0.5 + air*ka^0.5)^2
    end
end

function heatcapacity(params::SoilParams, totalWater, liquidWater, mineral, organic)
    @unpack _,_,_,_,_,cw,co,cm,ca,ci = params
    let air = 1.0 - totalWater - mineral - organic,
        ice = totalWater - liquidWater,
        water = liquidWater;
        water*cw + ice*ci + mineral*cm + organic*co + air*ca
    end
end

export thermalconductivity, heatcapacity
