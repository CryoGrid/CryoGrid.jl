const GridValues{Q,A} = NamedTuple{(:edges,:cells),NTuple{2,SubArray{Q,1,A,Tuple{StepRange{Int64,Int64}},true}}} where {Q,A<:AbstractArray}

"""
    struct Grid{S,G,Q,A} <: AbstractDiscretization{Q,1}

Represents the 1D spatial discretization on which time integration is performed. `S` is a `GridSpec`,
either `Edges` or `Cells` (always edges upon initial construction). The grid representation can be
converted (allocation free) between grid edges and cells via the `cells` and `edges` methods. `G`
represents the geometry/volume on which the vertical 1D discretization is applied. `A` is the underlying
array type, and `Q` is the numerical type (e.g. `Float64` or a `Unitful.Quantity`).
"""
struct Grid{S,G,Q,A} <: AbstractDiscretization{Q,1}
    geometry::G
    values::GridValues{Q,A}
    deltas::GridValues{Q,A}
    Grid(::Type{S}, values::GridValues{Q,A}, deltas::GridValues{Q,A}, geom::G) where {S<:GridSpec,Q,A,G<:Geometry} = new{S,G,Q,A}(geom,values,deltas)
    function Grid(vals::AbstractVector{Q}, geometry::G=UnitVolume()) where {G<:Geometry,Q<:Number}
        @assert issorted(vals) "grid values should be in ascending order"
        nedges = length(vals)
        ncells = nedges - 1 
        edges = copyto!(similar(vals), vals)
        cells = similar(edges, ncells)
        @. cells = (edges[1:end-1] + edges[2:end]) / (2*one(Q))
        # note that the distinction between edges and cells here is kind of redundant since the distance between
        # midpoints should always be the same as distance between boundaries...
        Δedges = edges[2:end] .- edges[1:end-1]
        Δcells = cells[2:end] .- cells[1:end-1]
        new{Edges,G,Q,typeof(edges)}(geometry,
            (edges=view(edges,1:1:nedges),cells=view(cells,1:1:ncells)),
            (edges=view(Δedges,1:1:nedges-1),cells=view(Δcells,1:1:ncells-1)),
        )
    end
    function Grid(grid::Grid{Edges,G,Q,A}, interval::ClosedInterval{Int}) where {G,Q,A}
        start, stop = interval.left, interval.right
        edges = @view grid.values.edges[start:stop]
        cells = @view grid.values.cells[start:stop-1]
        Δedges = @view grid.deltas.edges[start:stop-1]
        Δcells = @view grid.deltas.cells[start:stop-2]
        new{Edges,G,Q,A}(grid.geometry, (edges=edges,cells=cells), (edges=Δedges,cells=Δcells))
    end
    function Grid(grid::Grid{Cells,G,Q,A}, interval::ClosedInterval{Int}) where {G,Q,A}
        start, stop = interval.left, interval.right
        edges = @view grid.values.edges[start:stop+1]
        cells = @view grid.values.cells[start:stop]
        Δedges = @view grid.deltas.edges[start:stop]
        Δcells = @view grid.deltas.cells[start:stop-1]
        new{Cells,G,Q,A}(grid.geometry, (edges=edges,cells=cells), (edges=Δedges,cells=Δcells))
    end
    Grid(::Type{Cells}, grid::Grid{Edges,G,Q,A}) where {G,Q,A} = new{Cells,G,Q,A}(grid.geometry,grid.values,grid.deltas)
    Grid(::Type{Edges}, grid::Grid{Cells,G,Q,A}) where {G,Q,A} = new{Edges,G,Q,A}(grid.geometry,grid.values,grid.deltas)
end
ConstructionBase.constructorof(::Type{Grid{S,G,Q,A}}) where {S,G,Q,A} = (geom,values,deltas) -> Grid(S,values,deltas,geom)
Base.show(io::IO, grid::Grid{S,G}) where {S,G} = print(io, "Grid{$S}($(grid[1])..$(grid[end])) of length $(length(grid)) with geometry $G")
Base.show(io::IO, ::MIME{Symbol("text/plain")}, grid::Grid) = show(io, grid)

function subgridinds(grid::Grid{S,G,Q,A}, interval::Interval{L,R,Q}) where {S,G,Q,A,L,R}
    @assert interval.left <= interval.right "Invalid interval: $interval"
    vals = values(grid)
    # Determine indices which lie in the given interval
    l_ind = searchsortedfirst(vals, interval.left)
    r_ind = searchsortedlast(vals, interval.right)
    return (L == :closed ? l_ind : l_ind + 1)..(R == :closed ? r_ind : r_ind - 1)
end
@inline Δ(grid::Grid{Edges}) = grid.deltas.edges
@inline Δ(grid::Grid{Cells}) = grid.deltas.cells
@inline cells(grid::Grid{Edges}) = Grid(Cells, grid)
@inline cells(grid::Grid{Cells}) = grid
@inline edges(grid::Grid{Edges}) = grid
@inline edges(grid::Grid{Cells}) = Grid(Edges, grid)
@inline Base.values(grid::Grid{Edges}) = grid.values.edges
@inline Base.values(grid::Grid{Cells}) = grid.values.cells
@inline Base.similar(grid::Grid{Edges}) = Grid(copy(grid.values.edges))
@inline Base.similar(grid::Grid{Cells}) = similar(edges(grid)) |> cells
@inline Base.size(grid::Grid) = size(values(grid))
@inline Base.length(grid::Grid) = length(values(grid))
@propagate_inbounds Base.getindex(grid::Grid, i::Int) = values(grid)[i]
@propagate_inbounds Base.getindex(grid::Grid{S,G,Q,A}, interval::Interval{L,R,Q}) where {S,G,Q,A,L,R} = Grid(grid, subgridinds(grid,interval))
@propagate_inbounds Base.getindex(grid::Grid, interval::Interval{L,R,Int}) where {L,R} = Grid(grid, interval)
Base.setindex!(grid::Grid, args...) = error("setindex! is not allowed for Grid types")

# unit volume
@inline volume(grid::Grid{Cells,UnitVolume,Q}) where Q = Δ(edges(grid)).*oneunit(Q)^2
@inline area(::Grid{Edges,UnitVolume,Q}) where Q = oneunit(Q)^2
