"""
    SimpleSoil{Tpor,Tsat,Torg,Thp,Twp} <: SoilParameterization

Represents a simple organic/mineral soil mixutre in terms of its characteristic fractions:
i.e. natural porosity, saturation, and organic solid fraction. This is the standard CryoGrid representation
of a discrete soil volume.
"""
Base.@kwdef struct SimpleSoil{Tpor,Tsat,Torg,Thp,Twp} <: SoilParameterization
    por::Tpor = 0.5 # natural porosity
    sat::Tsat = 1.0 # saturation
    org::Torg = 0.0 # organic fraction of solid; mineral fraction is 1-org
    heat::Thp = SoilThermalProperties(SimpleSoil)
    water::Twp = SoilHydraulicProperties(SimpleSoil, fieldcapacity=0.20)
end

# Helper functions for obtaining soil compositions from characteristic fractions.
soilcomponent(::Val{var}, para::SimpleSoil) where var = soilcomponent(Val{var}(), para.por, para.sat, para.org)
soilcomponent(::Val{:θwi}, ϕ, θ, ω) = ϕ*θ
soilcomponent(::Val{:θp}, ϕ, θ, ω) = ϕ
soilcomponent(::Val{:θa}, ϕ, θ, ω) = ϕ*(1-θ)
soilcomponent(::Val{:θm}, ϕ, θ, ω) = (1-ϕ)*(1-ω)
soilcomponent(::Val{:θo}, ϕ, θ, ω) = (1-ϕ)*ω

# Soil methods (homogeneous)
saturation(soil::Soil{<:SimpleSoil}) = soil.para.sat
porosity(soil::Soil{<:SimpleSoil}) = soil.para.por
mineral(soil::Soil{<:SimpleSoil}) = soilcomponent(Val{:θm}(), soil.para)
organic(soil::Soil{<:SimpleSoil}) = soilcomponent(Val{:θo}(), soil.para)

# Soil thermal properties
SoilThermalProperties(
    ::Type{SimpleSoil};
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

# CryoGrid methods
CryoGrid.parameterize(para::SimpleSoil) = SimpleSoil(
    por = CryoGrid.parameterize(para.por, domain=0..1, desc="Natural porosity of the soil volume."),
    sat = CryoGrid.parameterize(para.sat, domain=0..1, desc="Initial water+ice saturation level of the soil volume."),
    org = CryoGrid.parameterize(para.org, domain=0..1, desc="Organic solid fraction of the soil volume."),
    heat = CryoGrid.parameterize(para.heat, domain=0..Inf),
    water = CryoGrid.parameterize(para.water, domain=0..Inf),
)

CryoGrid.variables(soil::Soil{<:SimpleSoil}) = CryoGrid.variables(soil, processes(soil))

CryoGrid.initializers(soil::Soil{<:SimpleSoil,THeat,<:WaterBalance}) where {THeat} = (
    initializer(:sat, soil.para.sat),
    initializers(soil, processes(soil))...,
)

# Heat

function Heat.thermalconductivities(soil::Soil{<:SimpleSoil})
    @unpack kh_w, kh_i, kh_a, kh_m, kh_o = thermalproperties(soil)
    return kh_w, kh_i, kh_a, kh_m, kh_o
end

function Heat.heatcapacities(soil::Soil{<:SimpleSoil})
    @unpack ch_w, ch_i, ch_a, ch_m, ch_o = thermalproperties(soil)
    return ch_w, ch_i, ch_a, ch_m, ch_o
end

"""
Gets the `ThermalProperties` for the given soil layer.
"""
Heat.thermalproperties(soil::Soil{<:SimpleSoil}) = soil.para.heat

# Hydrology

"""
Gets the `HydraulicProperties` for the given soil layer.
"""
Hydrology.hydraulicproperties(soil::Soil{<:SimpleSoil}) = soil.para.water

# field capacity
Hydrology.minwater(soil::Soil{<:SimpleSoil}, ::WaterBalance) = hydraulicproperties(soil).fieldcapacity

# porosity/max water
Hydrology.maxwater(soil::Soil{<:SimpleSoil}, water::WaterBalance) = porosity(soil, water)

# water content for soils without water balance
Hydrology.watercontent(soil::Soil{<:SimpleSoil,THeat,Nothing}, state) where {THeat} = soilcomponent(Val{:θwi}(), soil.para)
Hydrology.watercontent(soil::Soil{<:SimpleSoil,THeat,<:WaterBalance}, state) where {THeat} = state.θwi

# Heterogeneous

saturation(::Soil{<:Heterogeneous{<:SimpleSoil}}, state) = state.sat
porosity(::Soil{<:Heterogeneous{<:SimpleSoil}}, state) = state.por
mineral(::Soil{<:Heterogeneous{<:SimpleSoil}}, state) = state.θm
organic(::Soil{<:Heterogeneous{<:SimpleSoil}}, state) = state.θo

CryoGrid.variables(soil::Soil{<:Heterogeneous{<:SimpleSoil}}) = (
    # process variables
    variables(soil, processes(soil))...,
    # soil composition
    Diagnostic(:por, OnGrid(Cells), domain=0..1),
    Diagnostic(:sat, OnGrid(Cells), domain=0..1),
    Diagnostic(:org, OnGrid(Cells), domain=0..1),
    Diagnostic(:θwi, OnGrid(Cells), domain=0..1),
    Diagnostic(:θm, OnGrid(Cells), domain=0..1),
    Diagnostic(:θo, OnGrid(Cells), domain=0..1),
    # thermal properties
    Diagnostic(:kh_m, OnGrid(Cells), domain=0..Inf),
    Diagnostic(:kh_o, OnGrid(Cells), domain=0..Inf),
    Diagnostic(:kh_w, OnGrid(Cells), domain=0..Inf),
    Diagnostic(:kh_i, OnGrid(Cells), domain=0..Inf),
    Diagnostic(:kh_a, OnGrid(Cells), domain=0..Inf),
    Diagnostic(:ch_m, OnGrid(Cells), domain=0..Inf),
    Diagnostic(:ch_o, OnGrid(Cells), domain=0..Inf),
    Diagnostic(:ch_w, OnGrid(Cells), domain=0..Inf),
    Diagnostic(:ch_i, OnGrid(Cells), domain=0..Inf),
    Diagnostic(:ch_a, OnGrid(Cells), domain=0..Inf),
    # hydraulic properties
    Diagnostic(:kw_sat, OnGrid(Cells), domain=0..Inf),
    Diagnostic(:fieldcapacity, OnGrid(Cells), domain=0..1),
)

CryoGrid.initializers(soil::Soil{<:Heterogeneous{<:SimpleSoil}}) = (
    initializer(:por, map(para -> para.por, soil.para.profile)),
    initializer(:sat, map(para -> para.sat, soil.para.profile)),
    initializer(:org, map(para -> para.org, soil.para.profile)),
    initializer(:θwi, map(para -> para.por*para.sat, soil.para.profile)),
    initializer(:kh_m, map(para -> para.heat.kh_m, soil.para.profile)),
    initializer(:kh_o, map(para -> para.heat.kh_o, soil.para.profile)),
    initializer(:kh_w, map(para -> para.heat.kh_w, soil.para.profile)),
    initializer(:kh_i, map(para -> para.heat.kh_i, soil.para.profile)),
    initializer(:kh_a, map(para -> para.heat.kh_a, soil.para.profile)),
    initializer(:ch_m, map(para -> para.heat.ch_m, soil.para.profile)),
    initializer(:ch_o, map(para -> para.heat.ch_o, soil.para.profile)),
    initializer(:ch_w, map(para -> para.heat.ch_w, soil.para.profile)),
    initializer(:ch_i, map(para -> para.heat.ch_i, soil.para.profile)),
    initializer(:ch_a, map(para -> para.heat.ch_a, soil.para.profile)),
    initializer(:kw_sat, map(para -> para.water.kw_sat, soil.para.profile)),
    initializer(:fieldcapacity, map(para -> para.water.fieldcapacity, soil.para.profile)),
    initializers(soil, processes(soil))...,
)

function CryoGrid.initialcondition!(soil::Soil{<:Heterogeneous{<:SimpleSoil}}, state)
    # initialize θo and θm variables
    @. state.θo = soilcomponent(Val{:θo}(), state.por, state.sat, state.org)
    @. state.θm = soilcomponent(Val{:θm}(), state.por, state.sat, state.org)
    initialcondition!(soil, processes(soil), state)
end

Base.@propagate_inbounds function Heat.thermalproperties(::Soil{<:Heterogeneous{<:SimpleSoil}}, state, i)
    return ThermalProperties(
        kh_w = state.kh_w[i],
        kh_i = state.kh_i[i],
        kh_a = state.kh_a[i],
        kh_m = state.kh_m[i],
        kh_o = state.kh_o[i],
        ch_w = state.ch_w[i],
        ch_i = state.ch_i[i],
        ch_a = state.ch_a[i],
        ch_m = state.ch_m[i],
        ch_o = state.ch_o[i],
    )
end

Base.@propagate_inbounds function Heat.thermalconductivities(soil::Soil{<:Heterogeneous{<:SimpleSoil}}, state, i)
    @unpack kh_w, kh_i, kh_a, kh_m, kh_o = thermalproperties(soil, state, i)
    return (kh_w, kh_i, kh_a, kh_m, kh_o)
end

Base.@propagate_inbounds function Heat.heatcapacities(soil::Soil{<:Heterogeneous{<:SimpleSoil}}, state, i)
    @unpack ch_w, ch_i, ch_a, ch_m, ch_o = thermalproperties(soil, state, i)
    return (ch_w, ch_i, ch_a, ch_m, ch_o)
end

Base.@propagate_inbounds function Hydrology.hydraulicproperties(soil::Soil{<:Heterogeneous{<:SimpleSoil}}, state, i)
    return HydraulicProperties(
        kw_sat = state.kw_sat[i],
        fieldcapacity = state.fieldcapacity[i],
    )
end

Base.@propagate_inbounds Hydrology.minwater(soil::Soil{<:Heterogeneous{<:SimpleSoil}}, ::WaterBalance, state, i) = state.fieldcapacity[i]
