"""
    SnowpackParameterization

Base type for snowpack paramterization schemes.
"""
abstract type SnowpackParameterization <: CryoGrid.Parameterization end

"""
    Bulk{Tthresh} <: SnowpackParameterization

Simple, bulk snow scheme where snowpack is represented as a single grid cell with homogenous state.
"""
Base.@kwdef struct Bulk{Tthresh} <: SnowpackParameterization
    thresh::Tthresh = 0.05u"m" # snow threshold
end

"""
    Snowpack{Tpara<:SnowpackParameterization,Tsp} <: CryoGrid.SubSurface

Generic representation of a ground surface snow pack.
"""
Base.@kwdef struct Snowpack{Tpara<:SnowpackParameterization,Tsp} <: CryoGrid.SubSurface
    para::Tpara = Bulk()
    sp::Tsp = nothing
end
