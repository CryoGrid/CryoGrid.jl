# SURFEX paramterization
Base.@kwdef struct SURFEX{Tρs,Tρo,Tpo,Tsat,Twilt,Ttex<:SoilTexture,Tfc,Thp,Twp} <: SoilParameterization
    ρ_soc::Tρs = 65.0u"kg/m^3" # corresponds to 50% organic solid fraction
    ρ_org::Tρo = 1300u"kg/m^3" # pure organic matter density
    por_org::Tpo = 0.90 # organic material porosity
    sat::Tsat = 1.0 # saturation
    wilting_point::Twilt = 0.05 # hydraulic wilting point (minimum water holding capacity)
    texture::Ttex = SoilTexture(sand=0.50, clay=0.25) # soil texture
    freezecurve::Tfc = default_surfex_freezecurve(texture, ρ_soc, ρ_org, por_org)
    heat::Thp = SoilThermalProperties(SimpleSoil) # same thermal properties as SimpleSoil
    water::Twp = HydraulicProperties(fieldcapacity=0.20) # hydraulic properties
end

# SURFEX functions
surfex_organic(ρ_soc, ρ_org, por_org) = ρ_soc / ((1 - por_org)*ρ_org)

organic_porosity(surf::SURFEX) = surf.por_org

mineral_porosity(surf::SURFEX) = 0.49 - 0.11*surf.texture.sand

mineral_wilting_point(surf::SURFEX) = 0.37*sqrt(surf.texture.clay)

mineral_field_capacity(surf::SURFEX) = 0.45*surf.texture.clay^0.3496

organic_fraction(surf::SURFEX) = surfex_organic(surf.ρ_soc, surf.ρ_org, surf.por_org)

function default_surfex_freezecurve(texture::SoilTexture, ρ_soc, ρ_org, por_org)
    # TODO: find a more justified value for this threshold?
    # if soil is organic-rich, use organic freeze curve
    org = surfex_organic(ρ_soc, ρ_org, por_org)
    return if org > 0.10
        PainterKarra(swrc=VanGenuchten(:organic))
    else
        freezecurve(texture)
    end
end

# Soil methods
organic(soil::Soil{<:SURFEX}) = organic_fraction(soil.para)*(1 - porosity(soil))

mineral(soil::Soil{<:SURFEX}) = (1 - organic_fraction(soil.para))*(1 - porosity(soil))

saturation(soil::Soil{<:SURFEX}) = soil.para.sat

function porosity(soil::Soil{<:SURFEX})
    org = organic_fraction(soil.para)
    return (1-org)*mineral_porosity(soil.para) + org*organic_porosity(soil.para)
end

CryoGrid.initializers(soil::Soil{<:SURFEX,THeat,<:WaterBalance}) where {THeat} = (
    initializer(:sat, soil.para.sat),
    initializers(soil, processes(soil))...,
)

# Heat

Heat.freezecurve(soil::Soil{<:SURFEX}) = soil.para.freezecurve

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
