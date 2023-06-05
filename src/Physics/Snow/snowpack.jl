Base.@kwdef struct SnowpackProperties{Tmp,Thp,Twp}
    mass::Tmp = SnowMassProperties()
    heat::Thp = ThermalProperties()
    water::Twp = HydraulicProperties()
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
Base.@kwdef struct Snowpack{Tpara<:SnowpackParameterization,Tmass<:SnowMassBalance,Theat<:Optional{HeatBalance},Twater<:Optional{WaterBalance},Tprop,Taux} <: CryoGrid.SubSurface
    para::Tpara = Bulk()
    mass::Tmass = SnowMassBalance()
    heat::Theat = HeatBalance()
    water::Twater = nothing
    prop::Tprop = SnowpackProperties()
    aux::Taux = nothing
end

# type aliases for convenience
const PrescribedSnowpack{T} = Snowpack{T,<:SnowMassBalance{<:PrescribedSnow}} where {T}
const DynamicSnowpack{T} = Snowpack{T,<:SnowMassBalance{<:DynamicSnow}} where {T}

# Snow methods
"""
    threshold(::Snowpack)

Retrieves the snow cover threshold for this `Snowpack` layer to become active.
"""
threshold(::Snowpack) = 0.01 # meters
"""
    swe(::Snowpack, ::SnowMassBalance, state)

Retrieve the current snow water equivalent of the snowpack.
"""
swe(::Snowpack, ::SnowMassBalance, state) = state.swe
swe(::Snowpack, smb::SnowMassBalance{<:PrescribedSnow}, state) = smb.para.swe
swe(::Snowpack, smb::SnowMassBalance{<:PrescribedSnow{<:Forcing{u"m"}}}, state) = smb.para.swe(state.t)

"""
    snowdensity(::Snowpack, ::SnowMassBalance, state)

Retrieve the current snow density.
"""
snowdensity(::Snowpack, ::SnowMassBalance, state) = state.ρsn
snowdensity(::Snowpack, smb::SnowMassBalance{<:PrescribedSnow}, state) = smb.para.ρsn
snowdensity(::Snowpack, smb::SnowMassBalance{<:PrescribedSnow{Tswe,<:Forcing{u"kg/m^3"}}}, state) where {Tswe} = smb.para.ρsn(state.t)

"""
    snowdepth(::Snowpack, ::SnowMassBalance, state)

Retrieve the current snow depth.
"""
snowdepth(snow::Snowpack, ::SnowMassBalance, state) = CryoGrid.thickness(snow, state)

"""
    accumulation(snow::SnowMassBalance{<:DynamicSnow})

Get the snow accumulation scheme from the given `SnowMassBalance` parameterization.
"""
accumulation(snow::SnowMassBalance{<:DynamicSnow}) = snow.para.accumulation

"""
    ablation(snow::SnowMassBalance{<:DynamicSnow})

Get the snow ablation scheme from the given `SnowMassBalance` parameterization.
"""
ablation(snow::SnowMassBalance{<:DynamicSnow}) = snow.para.ablation

"""
    density(snow::SnowMassBalance{<:DynamicSnow})

Get the snow density scheme from the given `SnowMassBalance` parameterization.
"""
density(snow::SnowMassBalance{<:DynamicSnow}) = snow.para.density

snowvariables(::Snowpack) = (
    Diagnostic(:dsn, Scalar, u"m", domain=0..Inf),
    Diagnostic(:T_ub, Scalar, u"°C"),
)

# thermal properties of snowpack
Heat.thermalproperties(snow::Snowpack) = snow.prop.heat

# Default implementations of CryoGrid methods for Snowpack
CryoGrid.processes(snow::Snowpack{<:SnowpackParameterization,<:SnowMassBalance,<:HeatBalance,Nothing}) = Coupled(snow.mass, snow.heat)
CryoGrid.processes(snow::Snowpack{<:SnowpackParameterization,<:SnowMassBalance,<:HeatBalance,<:WaterBalance}) = Coupled(snow.mass, snow.water, snow.heat)

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

# interact! for coupled water/heat

function CryoGrid.interact!(top::Top, bc::WaterHeatBC, snow::Snowpack, ps::Coupled(SnowMassBalance, HeatBalance), stop, ssub)
    snowmass, heat = ps
    interact!(top, bc.water, snow, snowmass, stop, ssub)
    interact!(top, bc.heat, snow, heat, stop, ssub)
end

function CryoGrid.interact!(top::Top, bc::WaterHeatBC, snow::Snowpack, ps::Coupled(SnowMassBalance, WaterBalance, HeatBalance), stop, ssub)
    snowmass, water, heat = ps
    interact!(top, bc.water, snow, snowmass, stop, ssub)
    interact!(top, bc.water, snow, water, stop, ssub)
    interact!(top, bc.heat, snow, heat, stop, ssub)
end
