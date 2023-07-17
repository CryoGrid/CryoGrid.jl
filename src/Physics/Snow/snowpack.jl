# Snow mass balance
Utils.@properties SnowMassProperties(
    ρw = CryoGrid.Constants.ρw,
)

Base.@kwdef struct SnowpackProperties{Tmp,Thp,Twp}
    mass::Tmp = SnowMassProperties()
    heat::Thp = ThermalProperties()
    water::Twp = HydraulicProperties()
end

"""
    SnowAblationScheme

Base type for different snow ablation (i.e. melting or redistribution) schemes.
"""
abstract type SnowAblationScheme end

Base.@kwdef struct DegreeDayMelt{Tfactor,Tmax} <: SnowAblationScheme
    factor::Tfactor = 5.0u"mm/K/d"
    max_unfrozen::Tmax = 0.5
end

abstract type SnowAccumulationScheme end

Base.@kwdef struct LinearAccumulation{S} <: SnowAccumulationScheme
    rate_scale::S = 1.0 # scaling factor for snowfall rate
end

abstract type SnowDensityScheme end

# constant density (using Snowpack properties)
Base.@kwdef struct ConstantDensity{Tρsn}
    ρsn::Tρsn = 250.0u"kg/m^3"
end

abstract type SnowMassParameterization end

Base.@kwdef struct PrescribedSnow{Tswe,Tρsn} <: SnowMassParameterization
    swe::Tswe = 0.0u"m" # depth snow water equivalent [m]
    ρsn::Tρsn = 250.0u"kg/m^3" # snow density [kg m^-3]
end

Base.@kwdef struct DynamicSnow{TAcc,TAbl,TDen} <: SnowMassParameterization
    accumulation::TAcc = LinearAccumulation()
    ablation::TAbl = DegreeDayMelt()
    density::TDen = ConstantDensity()
end

Base.@kwdef struct SnowMassBalance{Tpara} <: CryoGrid.SubSurfaceProcess
    para::Tpara = DynamicSnow()
end

"""
    SnowpackParameterization

Base type for snowpack paramterization schemes.
"""
abstract type SnowpackParameterization <: CryoGrid.Parameterization end

"""
    Snowpack{Tpara<:SnowpackParameterization,Tprop,Taux} <: CryoGrid.SubSurface

Generic representation of a ground surface snow pack.
"""
Base.@kwdef struct Snowpack{Tpara<:SnowpackParameterization,Tmass<:SnowMassBalance,Theat<:Optional{HeatBalance},Twater<:WaterBalance,Tprop,Taux} <: CryoGrid.SubSurface
    para::Tpara = Bulk()
    mass::Tmass = SnowMassBalance()
    heat::Theat = HeatBalance()
    water::Twater = WaterBalance()
    prop::Tprop = SnowpackProperties()
    aux::Taux = nothing
end

# Snowpack parameterization type aliases
const PrescribedSnowpack{T} = Snowpack{T,<:SnowMassBalance{<:PrescribedSnow}} where {T}
const DynamicSnowpack{T} = Snowpack{T,<:SnowMassBalance{<:DynamicSnow}} where {T}

# Processes type aliases
const CoupledSnowWaterHeat{Tmass,Twater,Theat} = Coupled(SnowMassBalance, WaterBalance, HeatBalance)

# convenience type aliases
"""
    PrescribedSnowMassBalance{Tswe,Tρsn} = SnowMassBalance{PrescribedSnow{Tswe,Tρsn}} where {Tswe,Tρsn}
"""
const PrescribedSnowMassBalance{Tswe,Tρsn} = SnowMassBalance{PrescribedSnow{Tswe,Tρsn}} where {Tswe,Tρsn}
"""
    DynamicSnowMassBalance{TAcc,TAbl,TDen} = SnowMassBalance{DynamicSnow{TAcc,TAbl,TDen}} where {TAcc,TAbl,TDen}
"""
const DynamicSnowMassBalance{TAcc,TAbl,TDen} = SnowMassBalance{DynamicSnow{TAcc,TAbl,TDen}} where {TAcc,TAbl,TDen}

const SnowBC = BoundaryProcess{T} where {SnowMassBalance<:T<:SubSurfaceProcess}

# Heat methods;
# thermal properties of snowpack
Heat.thermalproperties(snow::Snowpack) = snow.prop.heat

# Hydrology methods;
Hydrology.hydraulicproperties(snow::Snowpack) = snow.prop.water

# max (fully saturated) water content
Hydrology.maxwater(::Snowpack, ::WaterBalance, state) = 1.0

# Default implementations of CryoGrid methods for Snowpack
CryoGrid.processes(snow::Snowpack) = Coupled(snow.mass, snow.water, snow.heat)

CryoGrid.thickness(::Snowpack, state, i::Integer=1) = abs(getscalar(state.Δz))

CryoGrid.midpoint(::Snowpack, state, i::Integer=1) = abs(getscalar(state.z) + getscalar(state.Δz)) / 2

CryoGrid.isactive(snow::Snowpack, state) = CryoGrid.thickness(snow, state) > threshold(snow)

CryoGrid.Volume(::Type{<:PrescribedSnowpack}) = CryoGrid.DiagnosticVolume()
CryoGrid.Volume(::Type{<:DynamicSnowpack}) = CryoGrid.PrognosticVolume()

# for prescribed snow depth/density, the mass balance is given so we do not need to do anything here
CryoGrid.computefluxes!(::Snowpack, ::SnowMassBalance{<:PrescribedSnow}, ssnow) = nothing

# volumetric fractions for snowpack
@inline function CryoGrid.volumetricfractions(::Snowpack, state, i)
    @inbounds let θwi = state.θwi[i],
        θw = state.θw[i],
        θa = 1.0 - θwi,
        θi = θwi - θw;
        return (θw, θi, θa)
    end
end

# always allow Top interactions with Snowpack
CryoGrid.caninteract(::Top, ::WaterBC, ::Snowpack, ::SnowMassBalance, s1, s2) = true

# default interact! for heat
function CryoGrid.interact!(
    top::Top,
    bc::HeatBC,
    snow::Snowpack,
    heat::HeatBalance,
    stop,
    ssnow
)
    @setscalar ssnow.T_ub = getscalar(stop.T_ub)
    # boundary flux
    ssnow.jH[1] += CryoGrid.boundaryflux(bc, top, heat, snow, stop, ssnow)
    return nothing
end
# default interact! for coupled water/heat
function CryoGrid.interact!(top::Top, bc::WaterHeatBC, snow::Snowpack, ps::Coupled(SnowMassBalance, HeatBalance), stop, ssub)
    waterbc, heatbc = bc
    snowmass, heat = ps
    interactmaybe!(top, waterbc, snow, snowmass, stop, ssub)
    interactmaybe!(top, heatbc, snow, heat, stop, ssub)
end
function CryoGrid.interact!(top::Top, bc::WaterHeatBC, snow::Snowpack, ps::Coupled(SnowMassBalance, WaterBalance, HeatBalance), stop, ssub)
    waterbc, heatbc = bc
    snowmass, water, heat = ps
    interactmaybe!(top, waterbc, snow, snowmass, stop, ssub)
    interactmaybe!(top, waterbc, snow, water, stop, ssub)
    interactmaybe!(top, heatbc, snow, heat, stop, ssub)
end
