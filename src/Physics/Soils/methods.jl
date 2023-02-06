# We use methods with optional index arguments `i` to allow for implementations both
# where these variables are treated as constants and as state variables.
# In the latter case, specializations should override only the index-free form
# and return a state vector instead of a scalar. The `getscalar` function will
# handle both the scalar and vector case!
"""
    mineral(soil::Soil, state, i)

Retrieves the mineral content for the given layer at grid cell `i`, if provided.
Defaults to using the scalar mineral content defined on `soil`.
"""
mineral(soil::Soil, state, i) = Utils.getscalar(mineral(soil, state), i)
mineral(::Soil, state) = state.θm

"""
    organic(soil::Soil, state, i)

Retrieves the organic content for the given layer at grid cell `i`, if provided.
Defaults to using the scalar organic content defined on `soil`.
"""
organic(soil::Soil, state, i) = Utils.getscalar(organic(soil, state), i)
organic(::Soil, state) = state.θo

"""
    porosity(soil::Soil, state, i)

Retrieves the porosity for the given layer at grid cell `i`, if provided.
Defaults to using the scalar porosity defined on `soil`.
"""
porosity(soil::Soil, state, i) = Utils.getscalar(porosity(soil, state), i)
porosity(::Soil, state) = state.θp

"""
    soilproperties(soil::Soil)

Retrieves the constant physical properties for this `Soil` layer. The default
implementation calls `soil.prop`. Note that this behavior can be overridden
as needed.
"""
soilproperties(soil::Soil) = soil.prop

"""
    soilproperties(para::SoilParameterization, proc::Process; prop_kwargs...)

Constructs a default set of soil properties/constants based on the given parameterization
and process(es). The default implementation simply constructs a `NamedTuple` from the given
keyword arguments.
"""
soilproperties(para::SoilParameterization, proc::Process; prop_kwargs...) = (; prop_kwargs...)
# Default behavior for coupled processes is to invoke SoilProperties on each individually and merge the results
soilproperties(para::SoilParameterization, procs::CoupledProcesses; prop_kwargs...) = reduce(merge, map(p -> SoilProperties(para, p; prop_kwargs...), procs))

"""
    SoilProfile(pairs::Pair{<:DistQuantity,<:SoilParameterization}...)

Alias for `Profile(pairs...)` assigning soil parameterizations to specific depths.
"""
SoilProfile(pairs::Pair{<:DistQuantity,<:SoilParameterization}...) = Profile(pairs...)
