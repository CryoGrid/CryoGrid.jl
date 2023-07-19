# SURFEX paramterization
Base.@kwdef struct SURFEX{Tρs,Tρo,Tpo,Ttex<:SoilTexture,Thp,Twp} <: SoilParameterization
    ρ_soc::Tρs = 65.0u"kg/m^3"
    ρ_org::Tρo = 1300u"kg/m^3"
    por_org::Tpo = 0.90
    wilting_point = 0.05
    texture::Ttex = SoilTexture()
    heat::Thp = SoilThermalProperties(MineralOrganic) # same thermal properties as MineralOrganic
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

function porosity(soil::Soil{<:SURFEX})
    org = organic(soil)
    por = org*soil.para.por_org + (1-org)*mineral_porosity(soil.para)
    return por
end

Heat.thermalproperties(soil::Soil{<:SURFEX}) = soil.para.heat

Hydrology.hydraulicproperties(soil::Soil{<:SURFEX}) = soil.para.water
