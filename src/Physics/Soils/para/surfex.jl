# SURFEX paramterization
Base.@kwdef struct SURFEX{Tρs,Tρo,Tpo,Tsat,Twilt,Ttex<:SoilTexture,Thp,Twp} <: SoilParameterization
    ρ_soc::Tρs = 65.0u"kg/m^3"
    ρ_org::Tρo = 1300u"kg/m^3"
    por_org::Tpo = 0.90
    sat::Tsat = 1.0
    wilting_point::Twilt = 0.05
    texture::Ttex = SoilTexture()
    heat::Thp = SoilThermalProperties(SimpleSoil) # same thermal properties as SimpleSoil
    water::Twp = HydraulicProperties(fieldcapacity=0.20) # hydraulic properties
end

# SURFEX functions
organic(surf::SURFEX) = surf.ρ_soc / ((1 - surf.por_org)*surf.ρ_org)

mineral_porosity(surf::SURFEX) = 0.49 - 0.11*surf.texture.sand

mineral_wilting_point(surf::SURFEX) = 0.37*sqrt(surf.texture.clay)

mineral_field_capacity(surf::SURFEX) = 0.45*surf.texture.clay^0.3496

# Soil methods
organic(soil::Soil{<:SURFEX}) = organic(soil.para)

mineral(soil::Soil{<:SURFEX}) = 1 - organic(soil.para)

saturation(soil::Soil{<:SURFEX}) = soil.para.sat

function porosity(soil::Soil{<:SURFEX})
    org = organic(soil)
    por = org*soil.para.por_org + (1-org)*mineral_porosity(soil.para)
    return por
end

CryoGrid.initializers(soil::Soil{<:SURFEX,THeat,<:WaterBalance}) where {THeat} = (
    initializer(:sat, soil.para.sat),
    initializers(soil, processes(soil))...,
)

# Heat

function Heat.thermalconductivities(soil::Soil{<:SURFEX})
    @unpack kh_w, kh_i, kh_a, kh_m, kh_o = thermalproperties(soil)
    return kh_w, kh_i, kh_a, kh_m, kh_o
end

function Heat.heatcapacities(soil::Soil{<:SURFEX})
    @unpack ch_w, ch_i, ch_a, ch_m, ch_o = thermalproperties(soil)
    return ch_w, ch_i, ch_a, ch_m, ch_o
end

"""
Gets the `ThermalProperties` for the given soil layer.
"""
Heat.thermalproperties(soil::Soil{<:SURFEX}) = soil.para.heat

# Soil thermal properties
SoilThermalProperties(
    ::Type{SURFEX};
    kh_w = ThermalProperties().kh_w,
    kh_i = ThermalProperties().kh_i,
    kh_a = ThermalProperties().kh_a,
    kh_o=0.25u"W/m/K", # organic [Hillel (1982)]
    kh_m=3.8u"W/m/K", # mineral [Hillel (1982)]
    ch_w = ThermalProperties().ch_w,
    ch_i = ThermalProperties().ch_i,
    ch_a = ThermalProperties().ch_a,
    ch_o=2.5e6u"J/K/m^3", # heat capacity organic
    ch_m=2.0e6u"J/K/m^3", # heat capacity mineral
    kwargs...,
) = ThermalProperties(; kh_w, kh_i, kh_a, kh_m, kh_o, ch_w, ch_i, ch_a, ch_m, ch_o, kwargs...)

# Hydrology

"""
Gets the `HydraulicProperties` for the given soil layer.
"""
Hydrology.hydraulicproperties(soil::Soil{<:SURFEX}) = soil.para.water

# field capacity
Hydrology.minwater(soil::Soil{<:SURFEX}, ::WaterBalance) = hydraulicproperties(soil).fieldcapacity

# porosity/max water
Hydrology.maxwater(soil::Soil{<:SURFEX}, water::WaterBalance) = porosity(soil)

# water content for soils without water balance
Hydrology.watercontent(soil::Soil{<:SURFEX,THeat,Nothing}, state) where {THeat} = porosity(soil)*saturation(soil, state)
