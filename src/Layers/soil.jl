"""
Represents the general type of the soil. Sand, Silt, and Clay are provided by default.
"""
abstract type SoilType end
struct Sand <: SoilType end
struct Silt <: SoilType end
struct Clay <: SoilType end
"""
Base type for soil parameterizations.
"""
abstract type SoilParameterization <: AbstractParameterization end
struct BySoilProperties <: SoilParameterization end
"""
Represents the composition of the soil in terms of fractions: excess ice, natural porosity, saturation, and organic/(mineral + organic).
"""
@with_kw struct SoilProperties <: Params @deftype Float64
    χ; @assert ϕ >= 0.0 && ϕ <= 1.0 # natural porosity
    ϕ; @assert θ >= 0.0 && θ <= 1.0 # saturation
    θ; @assert ω >= 0.0 && ω <= 1.0 # organic fraction of solid; mineral fraction is 1-ω
    ω; @assert χ >= 0.0 && χ <= 1.0 # excess ice fraction
end
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
    cw = 4.2e6xu"J/(K*m^3)" #[J/m^3K] heat capacity water
    co = 2.5e6xu"J/(K*m^3)" #[J/m^3K]  heat capacity organic
    cm = 2e6xu"J/(K*m^3)" #[J/m^3K]  heat capacity mineral
    ca = 0.00125e6xu"J/(K*m^3)" #[J/m^3K]  heat capacity pore space
    ci = 1.9e6xu"J/(K*m^3)" #[J/m^3K]  heat capacity ice
end
"""
Parameter type for Soil layers, includes thermal conductivity and heat capacity
constants as well as type/composition.
"""
@with_kw struct SoilParams{TType<:SoilType,TPara<:Parameterization,S} <: Params
    tc::SoilTCParams = SoilTCParams()
    hc::SoilHCParams = SoilHCParams()
    type::TType = Sand()
    para::TPara = nothing
    sp::S = nothing # user-defined specialization
end
"""
Basic Soil layer.
"""
struct Soil{TType,TPara,S} <: SubSurface
    profile::DimArray
    params::SoilParams{TType,TPara,S}
    function Soil(profile::DimArray; kwargs...)
        params = SoilParams(;kwargs...)
        new{typeof(params.type),typeof(params.para),typeof(params.sp)}(profile, params)
    end
end
"""
Alias/constructor for soil profile.
"""
function SoilProfile(vals::Pair{<:DistQuantity,SoilProperties}...)
    points = [d => tuple(props...) for (d,props) in vals]
    Profile(points...;names=fieldnames(SoilProperties))
end

export Soil, SoilProperties, SoilProfile, SoilParams, SoilType, Sand, Silt, Clay, SoilParameterization, BySoilProperties

# Helper functions for obtaining soil component fractions from soil properties.
soilcomp(::Val{:θx}, χ, ϕ, θ, ω) = χ
soilcomp(::Val{:θp}, χ, ϕ, θ, ω) = @. (1-χ)*ϕ*θ
soilcomp(::Val{:θm}, χ, ϕ, θ, ω) = @. (1-χ)*(1-ϕ)*(1-ω)
soilcomp(::Val{:θo}, χ, ϕ, θ, ω) = @. (1-χ)*(1-ϕ)*ω

Base.show(io::IO, soil::Soil{T,P}) where {T,P} = print(io, "Soil{$T,$P}($(soil.params))")

variables(::Soil{T,Nothing}) where T = (
    Diagnostic(:θw, Float64, OnGrid(Cells)),
    Diagnostic(:θp, Float64, OnGrid(Cells)),
    Diagnostic(:θx, Float64, OnGrid(Cells)),
    Diagnostic(:θl, Float64, OnGrid(Cells)),
    Diagnostic(:θm, Float64, OnGrid(Cells)),
    Diagnostic(:θo, Float64, OnGrid(Cells)),
)

variables(soil::Soil{T,BySoilProperties}) where T = (
    Diagnostic(:θw, Float64, OnGrid(Cells)),
    Diagnostic(:θp, Float64, OnGrid(Cells)),
    Diagnostic(:θx, Float64, OnGrid(Cells)),
    Diagnostic(:θl, Float64, OnGrid(Cells)),
    Diagnostic(:θm, Float64, OnGrid(Cells)),
    Diagnostic(:θo, Float64, OnGrid(Cells)),
    Parameter(:χ, collect(soil.profile[var=:χ]), 0..1),
    Parameter(:ϕ, collect(soil.profile[var=:ϕ]), 0..1),
    Parameter(:θ, collect(soil.profile[var=:θ]), 0..1),
    Parameter(:ω, collect(soil.profile[var=:ω]), 0..1),
)

function initialcondition!(soil::Soil{T,P}, state) where {T,P}
    # Helper functions for initializing soil composition state based on parameterization mode.
    getproperties(::Val{var}, soil::Soil{T,BySoilProperties}, state) where {var,T} = soilcomp(Val{var}(), state.params.χ, state.params.ϕ, state.params.θ, state.params.ω)
    getproperties(::Val{var}, soil::Soil{T,Nothing}, state) where {var,T} = soilcomp(Val{var}(), soil.profile[var=:χ], soil.profile[var=:ϕ], soil.profile[var=:θ], soil.profile[var=:ω])
    depths = length(size(soil.profile)) > 1 ? dims(soil.profile, :depth).val : [refdims(soil.profile)[1].val]
    for var in [:θx,:θp,:θm,:θo]
        arr = DimArray(similar(state[var], Union{Missing,eltype(state[var])}), (Z(state.grids[var]u"m",)))
        arr .= missing
        arr_sub = @view arr[Z(Near(depths))]
        arr_sub .= getproperties(Val{var}(), soil, state)
        Utils.ffill!(arr)
        state[var] .= skipmissing(arr)
    end
    @. state.θw = state.θx + state.θp
end
