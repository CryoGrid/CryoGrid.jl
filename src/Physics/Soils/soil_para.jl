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
    water::Twp = SoilHydraulicProperties(MineralOrganic, fieldcapacity=0.20)
end

# Helper functions for obtaining soil compositions from characteristic fractions.
soilcomponent(::Val{var}, para::MineralOrganic) where var = soilcomponent(Val{var}(), para.por, para.sat, para.org)
soilcomponent(::Val{:θwi}, ϕ, θ, ω) = ϕ*θ
soilcomponent(::Val{:θp}, ϕ, θ, ω) = ϕ
soilcomponent(::Val{:θa}, ϕ, θ, ω) = ϕ*(1-θ)
soilcomponent(::Val{:θm}, ϕ, θ, ω) = (1-ϕ)*(1-ω)
soilcomponent(::Val{:θo}, ϕ, θ, ω) = (1-ϕ)*ω

# Soil methods (homogeneous)
saturation(soil::Soil{<:MineralOrganic}) = soil.para.sat
porosity(soil::Soil{<:MineralOrganic}) = soil.para.por
mineral(soil::Soil{<:MineralOrganic}) = soilcomponent(Val{:θm}(), soil.para)
organic(soil::Soil{<:MineralOrganic}) = soilcomponent(Val{:θo}(), soil.para)

# Soil thermal properties
SoilThermalProperties(
    ::Type{MineralOrganic};
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
) = ThermalProperties(; kh_w, kh_i, kh_a, kh_m, kh_o, ch_w, ch_i, ch_a, ch_m, ch_o)

# CryoGrid methods
CryoGrid.parameterize(para::MineralOrganic) = MineralOrganic(
    por = CryoGrid.parameterize(para.por, domain=0..1, desc="Natural porosity of the soil volume."),
    sat = CryoGrid.parameterize(para.sat, domain=0..1, desc="Initial water+ice saturation level of the soil volume."),
    org = CryoGrid.parameterize(para.org, domain=0..1, desc="Organic solid fraction of the soil volume."),
    heat = CryoGrid.parameterize(para.heat, domain=0..Inf),
    water = CryoGrid.parameterize(para.water, domain=0..Inf),
)

CryoGrid.variables(soil::Soil{<:MineralOrganic}) = CryoGrid.variables(soil, processes(soil))

CryoGrid.initializers(soil::Soil{<:MineralOrganic,THeat,<:WaterBalance}) where {THeat} = (
    initializer(:sat, soil.para.sat),
    initializers(soil, processes(soil))...,
)

function Heat.thermalconductivities(soil::Soil{<:MineralOrganic})
    @unpack kh_w, kh_i, kh_a, kh_m, kh_o = thermalproperties(soil)
    return kh_w, kh_i, kh_a, kh_m, kh_o
end

function Heat.heatcapacities(soil::Soil{<:MineralOrganic})
    @unpack ch_w, ch_i, ch_a, ch_m, ch_o = thermalproperties(soil)
    return ch_w, ch_i, ch_a, ch_m, ch_o
end

"""
Gets the `ThermalProperties` for the given soil layer.
"""
Heat.thermalproperties(soil::Soil{<:MineralOrganic}) = soil.para.heat

"""
Gets the `HydraulicProperties` for the given soil layer.
"""
Hydrology.hydraulicproperties(soil::Soil{<:MineralOrganic}) = soil.para.water

# field capacity
Hydrology.minwater(soil::Soil{<:MineralOrganic}, ::WaterBalance) = hydraulicproperties(soil).fieldcapacity

# water content for soils without water balance
Hydrology.watercontent(soil::Soil{<:MineralOrganic,THeat,Nothing}, state) where {THeat} = soilcomponent(Val{:θwi}(), soil.para)
Hydrology.watercontent(soil::Soil{<:MineralOrganic,THeat,<:WaterBalance}, state) where {THeat} = state.θwi

"""
    Heterogeneous{V,N,D,Taux} <: SoilParameterization

Special `SoilParameterization` which wraps a `Profile` of another soil parameterization type
to indicate that it should be heterogeneous with over depth.
"""
Base.@kwdef struct Heterogeneous{V,N,D,Taux} <: SoilParameterization
    profile::SoilProfile{N,V,D}
    aux::Taux = nothing
    Heterogeneous(profile::SoilProfile{N,V,D}, aux::Taux=nothing) where {V,N,D,Taux} = new{V,N,D,Taux}(profile, aux)
end

Ground(soilprofile::SoilProfile; kwargs...) = Ground(Heterogeneous(soilprofile); kwargs...)

saturation(::Soil{<:Heterogeneous{<:MineralOrganic}}, state) = state.sat
porosity(::Soil{<:Heterogeneous{<:MineralOrganic}}, state) = state.por
mineral(::Soil{<:Heterogeneous{<:MineralOrganic}}, state) = state.θm
organic(::Soil{<:Heterogeneous{<:MineralOrganic}}, state) = state.θo

CryoGrid.variables(soil::Soil{<:Heterogeneous{<:MineralOrganic}}) = (
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

CryoGrid.initializers(soil::Soil{<:Heterogeneous{<:MineralOrganic}}) = (
    initializer(:por, map(para -> para.por, soil.para.profile)),
    initializer(:sat, map(para -> para.sat, soil.para.profile)),
    initializer(:org, map(para -> para.org, soil.para.profile)),
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

function CryoGrid.initialcondition!(soil::Soil{<:Heterogeneous{<:MineralOrganic}}, state)
    # initialize θo and θm variables
    @. state.θo = soilcomponent(Val{:θo}(), state.por, state.sat, state.org)
    @. state.θm = soilcomponent(Val{:θm}(), state.por, state.sat, state.org)
    initialcondition!(soil, processes(soil), state)
end

Base.@propagate_inbounds function Heat.thermalproperties(::Soil{<:Heterogeneous{<:MineralOrganic}}, state, i)
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

Base.@propagate_inbounds function Heat.thermalconductivities(soil::Soil{<:Heterogeneous{<:MineralOrganic}}, state, i)
    @unpack kh_w, kh_i, kh_a, kh_m, kh_o = thermalproperties(soil, state, i)
    return (kh_w, kh_i, kh_a, kh_m, kh_o)
end

Base.@propagate_inbounds function Heat.heatcapacities(soil::Soil{<:Heterogeneous{<:MineralOrganic}}, state, i)
    @unpack ch_w, ch_i, ch_a, ch_m, ch_o = thermalproperties(soil, state, i)
    return (ch_w, ch_i, ch_a, ch_m, ch_o)
end

Base.@propagate_inbounds function Hydrology.hydraulicproperties(soil::Soil{<:Heterogeneous{<:MineralOrganic}}, state, i)
    return HydraulicProperties(
        kw_sat = state.kw_sat[i],
        fieldcapacity = state.fieldcapacity[i],
    )
end

Base.@propagate_inbounds Hydrology.minwater(soil::Soil{<:Heterogeneous{<:MineralOrganic}}, ::WaterBalance, state, i) = state.fieldcapacity[i]
