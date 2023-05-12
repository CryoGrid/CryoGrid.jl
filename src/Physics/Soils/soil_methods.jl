# Soils module methods
# We use methods with optional index arguments `i` to allow for implementations both
# where these variables are treated as constants and as state variables.
# In the latter case, specializations should override only the index-free form
# and return a state vector instead of a scalar. The `getscalar` function will
# handle both the scalar and vector case!
"""
    mineral(soil::Soil, state, i)
    mineral(soil::Soil, state)
    mineral(soil::Soil)

Retrieves the mineral content for the given layer at grid cell `i`, if provided.
"""
mineral(soil::Soil, state, i) = Utils.getscalar(mineral(soil, state), i)
mineral(soil::Soil, state) = mineral(soil)
mineral(soil::Soil) = error("not implemented for $(typeof(soil))")

"""
    organic(soil::Soil, state, i)
    organic(soil::Soil, state)
    organic(soil::Soil)

Retrieves the organic content for the given layer at grid cell `i`, if provided.
"""
organic(soil::Soil, state, i) = Utils.getscalar(organic(soil, state), i)
organic(soil::Soil, state) = mineral(soil)
organic(::Soil) = error("not implemented for $(typeof(soil))")

"""
    porosity(soil::Soil, state, i)
    porosity(soil::Soil, state)
    porosity(soil::Soil)

Retrieves the porosity for the given layer at grid cell `i`, if provided.
"""
porosity(soil::Soil, state, i) = Utils.getscalar(porosity(soil, state), i)
porosity(soil::Soil, state) = porosity(soil)
porosity(::Soil) = error("not implemented for $(typeof(soil))")

"""
    saturation(soil, state, i)
    saturation(soil::Soil, state)
    saturation(soil::Soil)

Retrieves the saturation level for the given layer at grid cell `i`, if provided.
"""
saturation(soil::Soil, state, i) = Utils.getscalar(saturation(soil, state), i)
saturation(soil::Soil, state) = saturation(soil)
# For soil layers with water balance, use saturation state variable
saturation(::Soil{<:Any,<:Any,<:WaterBalance}, state) = state.sat
saturation(::Soil) = error("not implemented for $(typeof(soil))")

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
