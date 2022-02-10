const GridValues{A} = NamedTuple{(:edges,:cells),NTuple{2,A}} where {A<:AbstractVector}

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
    values::GridValues{A}
    deltas::GridValues{A}
    bounds::UnitRange{Int}
    Grid(::Type{S}, values::GridValues{A}, deltas::GridValues{A}, geom::G, bounds::UnitRange{Int}=1:length(values)) where {S<:GridSpec,Q,A<:AbstractVector{Q},G<:Geometry} = new{S,G,Q,A}(geom,values,deltas,bounds)
    function Grid(vals::AbstractVector{Q}, geometry::G=UnitVolume()) where {G<:Geometry,Q<:Number}
        @assert issorted(vals) "grid values should be in ascending order"
        nedges = length(vals)
        ncells = nedges - 1
        edges = copyto!(similar(vals), vals)
        cells = similar(edges, ncells)
        @. cells = (edges[1:end-1] + edges[2:end]) / (2*one(Q))
        Δedges = edges[2:end] .- edges[1:end-1]
        Δcells = cells[2:end] .- cells[1:end-1]
        new{Edges,G,Q,typeof(edges)}(geometry, (edges=edges, cells=cells), (edges=Δedges, cells=Δcells), 1:length(edges))
    end
    function Grid(grid::Grid{S,G,Q,A}, interval::ClosedInterval{Int}) where {S,G,Q,A}
        new{S,G,Q,A}(grid.geometry, grid.values, grid.deltas, leftendpoint(interval):rightendpoint(interval))
    end
    Grid(::Type{Cells}, grid::Grid{Edges,G,Q,A}) where {G,Q,A} = new{Cells,G,Q,A}(grid.geometry, grid.values, grid.deltas, grid.bounds)
    Grid(::Type{Edges}, grid::Grid{Cells,G,Q,A}) where {G,Q,A} = new{Edges,G,Q,A}(grid.geometry, grid.values, grid.deltas, grid.bounds)
end
ConstructionBase.constructorof(::Type{Grid{S,G,Q,A}}) where {S,G,Q,A} = (geom,values,deltas,bounds) -> Grid(S,values,deltas,geom,bounds)
Base.show(io::IO, ::MIME"text/plain", grid::Grid) = show(io, grid)
function Base.show(io::IO, grid::Grid{S,G}) where {S,G}
    if length(grid) == length(grid.values.edges)
        print(io, "Grid{$S}($(grid[1])..$(grid[end])) of length $(length(grid)) with geometry $G")
    else
        print(io, "Grid{$S}($(grid[1])..$(grid[end])) of length $(length(grid)) (child of Grid{$S}$(parent(grid)[1])..$(parent(grid)[end]) of length $(length(parent(grid)))) with geometry $G")
    end
end

function subgridinds(grid::Grid, interval::Interval{L,R}) where {L,R}
    @assert interval.left <= interval.right "Invalid interval: $interval"
    # Determine indices which lie in the given interval
    l_ind = searchsortedfirst(grid, interval.left)
    r_ind = searchsortedlast(grid, interval.right)
    return (L == :closed ? l_ind : l_ind + 1)..(R == :closed ? r_ind : r_ind - 1)
end
@inline bounds(grid::Grid{Edges}) = grid.bounds
@inline bounds(grid::Grid{Cells}) = first(grid.bounds):last(grid.bounds)-1
@inline Δbounds(grid::Grid{Edges}) = first(grid.bounds):last(grid.bounds)-1
@inline Δbounds(grid::Grid{Cells}) = first(grid.bounds):last(grid.bounds)-2
@inline Δ(grid::Grid{Edges}) = view(grid.deltas.edges, Δbounds(grid))
@inline Δ(grid::Grid{Cells}) = view(grid.deltas.cells, Δbounds(grid))
@inline cells(grid::Grid{Edges}) = Grid(Cells, grid)
@inline cells(grid::Grid{Cells}) = grid
@inline edges(grid::Grid{Edges}) = grid
@inline edges(grid::Grid{Cells}) = Grid(Edges, grid)
@inline Base.parent(grid::Grid{Edges}) = Grid(grid, 1..length(grid.values.edges))
@inline Base.parent(grid::Grid{Cells}) = Grid(grid, 1..length(grid.values.cells))
@inline Base.collect(grid::Grid) = collect(values(grid))
@inline Base.values(grid::Grid{Edges}) = view(grid.values.edges, bounds(grid))
@inline Base.values(grid::Grid{Cells}) = view(grid.values.cells, bounds(grid))
@inline Base.similar(grid::Grid{Edges}) = Grid(copy(grid.values.edges))
@inline Base.similar(grid::Grid{Cells}) = similar(edges(grid)) |> cells
@inline Base.size(grid::Grid) = (length(grid),)
@inline Base.length(grid::Grid) = last(bounds(grid)) - first(bounds(grid)) + 1
@inline Base.firstindex(grid::Grid) = 1
@inline Base.lastindex(grid::Grid) = length(grid)
@propagate_inbounds Base.getindex(grid::Grid, i::Int) = values(grid)[i]
@propagate_inbounds Base.getindex(grid::Grid, i::AbstractRange) = grid[grid[first(i)]..grid[last(i)]]
@propagate_inbounds Base.getindex(grid::Grid{S,G,Q,A}, interval::Interval{L,R,Q}) where {S,G,Q,A,L,R} = grid[subgridinds(grid, interval)]
@propagate_inbounds Base.getindex(grid::Grid, interval::Interval{L,R,Int}) where {L,R} = Grid(grid, first(grid.bounds)+interval.left-1..first(grid.bounds)+interval.right-1)
Base.setindex!(::Grid, args...) = error("setindex! is not allowed for Grid types")
"""
    updategrid!(grid::Grid{Edges,G,Q}, vals::Q) where {G,Q}

Overwrites `grid` edges with `vals`, and recomputes grid centers/deltas to be consistent with the new grid.
"""
function updategrid!(grid::Grid{Edges,G,Q}, vals::AbstractVector{Q}) where {G,Q}
    z_edges = values(grid)
    z_cells = values(cells(grid))
    Δz_edges = Δ(grid)
    Δz_cells = Δ(cells(grid))
    z_edges .= vals
    z_cells .= (z_edges[1:end-1] .+ z_edges[2:end]) ./ (2*one(Q))
    Δz_edges .= z_edges[2:end] .- z_edges[1:end-1]
    Δz_cells .= z_cells[2:end] .- z_cells[1:end-1]
    @assert issorted(parent(grid)) "updated grid values are invalid; grid edges must be strictly non-decreasing"
    return grid
end

# unit volume
@inline volume(grid::Grid{Cells,UnitVolume,Q}) where Q = Δ(edges(grid)).*oneunit(Q)^2
@inline area(::Grid{Edges,UnitVolume,Q}) where Q = oneunit(Q)^2
