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

# homogeneous
saturation(soil::Soil{<:MineralOrganic}) = soil.para.sat
porosity(soil::Soil{<:MineralOrganic}) = soil.para.por
mineral(soil::Soil{<:MineralOrganic}) = soilcomponent(Val{:θm}(), soil.para)
organic(soil::Soil{<:MineralOrganic}) = soilcomponent(Val{:θo}(), soil.para)

# heterogeneous (soil composition discretized on the grid)
saturation(::Soil{<:Heterogeneous{<:MineralOrganic}}, state) = state.sat
porosity(::Soil{<:Heterogeneous{<:MineralOrganic}}, state) = state.por
mineral(::Soil{<:Heterogeneous{<:MineralOrganic}}, state) = state.θm
organic(::Soil{<:Heterogeneous{<:MineralOrganic}}, state) = state.θo

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

# CryoGrid methods
CryoGrid.parameterize(para::MineralOrganic) = MineralOrganic(
    por = CryoGrid.parameterize(para.por, domain=0..1),
    sat = CryoGrid.parameterize(para.sat, domain=0..1),
    org = CryoGrid.parameterize(para.org, domain=0..1),
)

CryoGrid.variables(soil::Soil{<:MineralOrganic}) = CryoGrid.variables(soil, processes(soil))
CryoGrid.variables(soil::Soil{<:Heterogeneous{<:MineralOrganic}}) = (
    Diagnostic(:por, OnGrid(Cells), domain=0..1),
    Diagnostic(:sat, OnGrid(Cells), domain=0..1),
    Diagnostic(:org, OnGrid(Cells), domain=0..1),
    Diagnostic(:θwi, OnGrid(Cells), domain=0..1),
    Diagnostic(:θm, OnGrid(Cells), domain=0..1),
    Diagnostic(:θo, OnGrid(Cells), domain=0..1),
    variables(soil, processes(soil))...
)

CryoGrid.initializers(soil::Soil{<:MineralOrganic,THeat,<:WaterBalance}) where {THeat} = (
    initializer(:sat, soil.para.sat),
    initializers(soil, processes(soil))...,
)

CryoGrid.initializers(soil::Soil{<:Heterogeneous{<:MineralOrganic}}) = (
    initializer(:por, soil.para.para.por),
    initializer(:sat, soil.para.para.sat),
    initializer(:org, soil.para.para.org),
    initializer(:θwi, (soil,state) -> state.θwi .= state.por.*state.sat),
    initializers(soil, processes(soil))...,
)

function CryoGrid.initialcondition!(soil::Soil{<:Heterogeneous{<:MineralOrganic}}, state)
    # initialize θo and θm variables
    @. state.θo = soilcomponent(Val{:θo}(), state.por, state.sat, state.org)
    @. state.θm = soilcomponent(Val{:θm}(), state.por, state.sat, state.org)
    initialcondition!(soil, processes(soil), state)
end

"""
Gets the `ThermalProperties` for the given soil layer.
"""
Heat.thermalproperties(soil::Soil{<:MineralOrganic}) = soil.para.heat
Heat.thermalproperties(soil::Soil{<:Heterogeneous{<:MineralOrganic}}) = soil.para.para.heat

"""
Gets the `HydraulicProperties` for the given soil layer.
"""
Hydrology.hydraulicproperties(soil::Soil{<:MineralOrganic}) = soil.para.water
Hydrology.hydraulicproperties(soil::Soil{<:Heterogeneous{<:MineralOrganic}}) = soil.para.para.water

# water content for soils without water balance
Hydrology.watercontent(soil::Soil{<:MineralOrganic,THeat,Nothing}, state) where {THeat} = soilcomponent(Val{:θwi}(), soil.para)
Hydrology.watercontent(soil::Soil{<:MineralOrganic,THeat,<:WaterBalance}, state) where {THeat} = state.θwi