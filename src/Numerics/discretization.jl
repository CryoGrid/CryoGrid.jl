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

abstract type AutoGridSpacing end
Base.@kwdef struct LinearSpacing <: AutoGridSpacing
    min_thick = 0.02u"m" # minimum cell thickness
    max_cells_per_layer = 100
end

"""
    AutoGrid <: DiscretizationStrategy
"""
Base.@kwdef struct AutoGrid{S<:AutoGridSpacing} <: DiscretizationStrategy
    spacing::S = LinearSpacing()
    z0 = 0.0u"m" # minimum/initial depth
    geometry = UnitRectangle() # lateral geometry of the volume 
end

# grid initialization
"""
    makegrid(strategy::DiscretizationStrategy, bounds::NTuple{2,<:DistQuantity})
    makegrid(layer::Layer, strategy::DiscretizationStrategy, bounds::NTuple{2,<:DistQuantity})

Constructs a `Grid` spanning `bounds` using the given `strategy`. `makegrid` can also be specified
for specific `Layer` types when layers have specific discretization requirements.
"""
function makegrid(strategy::PresetGrid, bounds::NTuple{2,<:DistQuantity})
    return strategy.grid[bounds[1]..bounds[2]]
end
function makegrid(strategy::AutoGrid{LinearSpacing}, bounds::NTuple{2,T}) where {T<:DistQuantity}
    @unpack spacing, z0 = strategy
    @unpack min_thick, max_cells_per_layer = spacing
    z1, z2 = bounds
    dz = min_thick
    n = min(max_cells_per_layer, Int(ceil((z2 - z1) / dz)))
    gridvals = LinRange(z1, z2, n)
    return Grid(round.(T, gridvals, digits=12))
end
makegrid(::Layer, strategy, bounds) = makegrid(strategy, bounds)

# state variable instantiation
"""
    instantiate(::Var, ::T, ::Type{A}) where {T,N,D<:AbstractDiscretization{T,N},A<:AbstractArray{T,N}}

Produces a discretization of the given variable based on `T` and array type `A`.
"""
instantiate(var::Var, d::AbstractDiscretization{Q,N}) where {Q,N} = instantiate(var, d, Array{vartype(var),N})
instantiate(var::Var, grid::Grid, ::Type{A}) where {A<:AbstractVector} = zero(similar(A{vartype(var)}, dimlength(var.dim, length(edges(grid)))))
