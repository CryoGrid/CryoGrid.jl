abstract type AbstractDiscretization{Q,N} <: DenseArray{Q,N} end

abstract type Geometry end
struct UnitRectangle <: Geometry end

"""
Base type for discretization "strategies" that generate a spatial discretization (typically a grid).
"""
abstract type DiscretizationStrategy end
"""
Simple discretization strategy that just supplies a pre-specified grid.
"""
struct PresetGrid{TGrid<:AbstractDiscretization} <: DiscretizationStrategy
    grid::TGrid
end
"""
    AutoGrid <: DiscretizationStrategy
"""
Base.@kwdef struct AutoGrid <: DiscretizationStrategy
    min_thick = 0.02u"m" # minimum cell thickness
    max_cells_per_layer = 100
    z0 = 0.0u"m" # minimum/initial depth
    geometry = UnitRectangle() # lateral geometry of the volume 
end
