# Note that here "soil_types" refers to programming types, not geological soil types!

"""
    SoilParameterization

Base type for parameterizations of soil consituents.
"""
abstract type SoilParameterization <: GroundParameterization end

"""
    Soil{Tpara,Theat,Twater}

Type alias for any `AbstractGround` layer with parameterization of type `SoilParameterization`.
"""
const Soil{Tpara,Theat,Twater} = AbstractGround{Tpara,Theat,Twater} where {Tpara<:SoilParameterization,Theat<:Optional{HeatBalance},Twater<:Optional{WaterBalance}}

"""
    SoilProfile{N,V,D} = Profile{N,V,D} where {N,V<:SoilParameterization,D<:DistQuantity}

Alias for depthwise `Profile` where the values are `SoilParameterization` types.
"""
const SoilProfile{N,V,D} = Profile{N,V,D} where {N,V<:SoilParameterization,D<:DistQuantity}

# Constructors
"""
    SoilProfile(pairs::Pair{<:DistQuantity,<:SoilParameterization}...)

Alias for `Profile(pairs...)` specific for `SoilProfile`s.
"""
SoilProfile(pairs::Pair{<:DistQuantity,<:SoilParameterization}...) = Profile(pairs...)
