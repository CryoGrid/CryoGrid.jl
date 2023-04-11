"""
    SoilParameterization

Base type for parameterizations of soil consituents.
"""
abstract type SoilParameterization end

"""
    Soil{Tpara,Theat<:Optional{HeatBalance},Twater<:Optional{WaterBalance}} <: SubSurface

Base type for all soil or soil-like layers.
"""
abstract type Soil{Tpara,Theat<:Optional{HeatBalance},Twater<:Optional{WaterBalance}} <: SubSurface end

Base.@kwdef struct SoilProperties{Thp,Twp}
    heat::Thp = SoilThermalProperties()
    water::Twp = HydraulicProperties()
end

"""
    HomogeneousSoil{Tpara<:SoilParameterization,Theat<:Optional{HeatBalance},Twater<:Optional{WaterBalance},Tsp,Tprop} <: Soil{Tpara,Theat,Twater}

Generic, homogeneous Soil layer, i.e. material is assumed to be uniformly mixed.
"""
Base.@kwdef struct HomogeneousSoil{Tpara<:SoilParameterization,Theat<:Optional{HeatBalance},Twater<:Optional{WaterBalance},Tsp,Tprop} <: Soil{Tpara,Theat,Twater}
    para::Tpara = HomogeneousMixture() # soil parameterization
    prop::Tprop = SoilProperties() # soil properties
    heat::Theat = HeatBalance() # heat conduction
    water::Twater = nothing # water balance
    sp::Tsp = nothing # user-defined specialization
end

# Soils module methods
# We use methods with optional index arguments `i` to allow for implementations both
# where these variables are treated as constants and as state variables.
# In the latter case, specializations should override only the index-free form
# and return a state vector instead of a scalar. The `getscalar` function will
# handle both the scalar and vector case!
"""
    mineral(soil::Soil, state, i)

Retrieves the mineral content for the given layer at grid cell `i`, if provided.
"""
mineral(soil::Soil, state, i) = Utils.getscalar(mineral(soil, state), i)
mineral(::Soil, state) = state.θm

"""
    organic(soil::Soil, state, i)

Retrieves the organic content for the given layer at grid cell `i`, if provided.
"""
organic(soil::Soil, state, i) = Utils.getscalar(organic(soil, state), i)
organic(::Soil, state) = state.θo

"""
    porosity(soil::Soil, state, i)

Retrieves the porosity for the given layer at grid cell `i`, if provided.
"""
porosity(soil::Soil, state, i) = Utils.getscalar(porosity(soil, state), i)
porosity(::Soil, state) = state.θsat

"""
    saturation(soil, state, i)

Retrieves the saturation level for the given layer at grid cell `i`, if provided.
"""
saturation(soil::Soil, state, i) = Utils.getscalar(saturation(soil, state), i)
saturation(::Soil, state) = state.sat

"""
    soilproperties(soil::Soil)

Retrieves the constant physical properties for this `Soil` layer. The default
implementation calls `soil.prop`. Note that this behavior can be overridden
as needed.
"""
soilproperties(soil::Soil) = soil.prop

# Constructors

"""
    SoilProfile(pairs::Pair{<:DistQuantity,<:SoilParameterization}...)

Alias for `Profile(pairs...)` assigning soil parameterizations to specific depths.
"""
SoilProfile(pairs::Pair{<:DistQuantity,<:SoilParameterization}...) = Profile(pairs...)

# CryoGrid core methods

CryoGrid.processes(soil::Soil{<:SoilParameterization,<:HeatBalance,Nothing}) = soil.heat
CryoGrid.processes(soil::Soil{<:SoilParameterization,Nothing,<:WaterBalance}) = soil.water
CryoGrid.processes(soil::Soil{<:SoilParameterization,<:HeatBalance,<:WaterBalance}) = Coupled(soil.water, soil.heat)
