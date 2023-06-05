"""
    MineralOrganic{Tpor,Tsat,Torg} <: SoilParameterization

Represents a simple organic/mineral soil mixutre in terms of its characteristic fractions:
i.e. natural porosity, saturation, and organic solid fraction. This is the standard CryoGrid representation
of a discrete soil volume.
"""
Base.@kwdef struct MineralOrganic{Tpor,Tsat,Torg,Thp,Twp} <: SoilParameterization
    por::Tpor = 0.5 # natural porosity
    sat::Tsat = 1.0 # saturation
    org::Torg = 0.0 # organic fraction of solid; mineral fraction is 1-org
    heat::Thp = SoilThermalProperties(MineralOrganic)
    water::Twp = HydraulicProperties()
end

# Helper functions for obtaining soil compositions from characteristic fractions.
soilcomponent(::Val{var}, para::MineralOrganic) where var = soilcomponent(Val{var}(), para.por, para.sat, para.org)
soilcomponent(::Val{:θwi}, ϕ, θ, ω) = ϕ*θ
soilcomponent(::Val{:θp}, ϕ, θ, ω) = ϕ
soilcomponent(::Val{:θa}, ϕ, θ, ω) = ϕ*(1-θ)
soilcomponent(::Val{:θm}, ϕ, θ, ω) = (1-ϕ)*(1-ω)
soilcomponent(::Val{:θo}, ϕ, θ, ω) = (1-ϕ)*ω

# Soil methods
saturation(soil::Soil{<:MineralOrganic}) = soil.para.sat
porosity(soil::Soil{<:MineralOrganic}) = soil.para.por
mineral(soil::Soil{<:MineralOrganic}) = soilcomponent(Val{:θm}(), soil.para)
organic(soil::Soil{<:MineralOrganic}) = soilcomponent(Val{:θo}(), soil.para)

# CryoGrid methods
CryoGrid.parameterize(para::MineralOrganic) = MineralOrganic(
    por = CryoGrid.parameterize(para.por, domain=0..1),
    sat = CryoGrid.parameterize(para.sat, domain=0..1),
    org = CryoGrid.parameterize(para.org, domain=0..1),
)
CryoGrid.variables(soil::Soil{<:MineralOrganic}) = CryoGrid.variables(soil, processes(soil))

# Soil thermal properties
const DefaultThermalProperties = Heat.ThermalProperties()
SoilThermalProperties(
    ::Type{MineralOrganic};
    kh_w = DefaultThermalProperties.kh_w,
    kh_i = DefaultThermalProperties.kh_i,
    kh_a = DefaultThermalProperties.kh_a,
    kh_o=0.25u"W/m/K", # organic [Hillel (1982)]
    kh_m=3.8u"W/m/K", # mineral [Hillel (1982)]
    ch_w = DefaultThermalProperties.ch_w,
    ch_i = DefaultThermalProperties.ch_i,
    ch_a = DefaultThermalProperties.ch_a,
    ch_o=2.5e6u"J/K/m^3", # heat capacity organic
    ch_m=2.0e6u"J/K/m^3", # heat capacity mineral
) = ThermalProperties(; kh_w, kh_i, kh_a, kh_m, kh_o, ch_w, ch_i, ch_a, ch_m, ch_o)

"""
Gets the `ThermalProperties` for the given soil layer.
"""
Heat.thermalproperties(soil::Soil{<:MineralOrganic}) = soil.para.heat

"""
Gets the `HydraulicProperties` for the given soil layer.
"""
Hydrology.hydraulicproperties(soil::Soil{<:MineralOrganic}) = soil.para.water

# water content for soils without water balance
Hydrology.watercontent(soil::Soil{<:MineralOrganic,THeat,Nothing}, state) where {THeat} = soilcomponent(Val{:θwi}(), soil.para)
Hydrology.watercontent(soil::Soil{<:MineralOrganic,THeat,<:WaterBalance}, state) where {THeat} = state.θwi
