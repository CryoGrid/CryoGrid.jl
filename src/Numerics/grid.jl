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
    l_ind = searchsortedlast(grid, interval.left)
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
Base.setindex!(grid::Grid{Edges}, val, i...) = setindex!(values(grid), val, i...)
Base.setindex!(::Grid{Cells}, args...) = error("setindex! is permitted only for edge grids; use `edges(grid)` and call `updategrid!` directly after.")
"""
    updategrid!(grid::Grid{Edges,G,Q}, vals::Q=grid) where {G,Q}

Overwrites `grid` edges with `vals`, and recomputes grid centers/deltas to be consistent with the new grid.
"""
function updategrid!(grid::Grid{Edges,G,Q}, vals::AbstractVector{Q}=grid) where {G,Q}
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

# grid discretizations
"""
    discretize([::Type{A}], ::T,  ::Var) where {T,N,D<:AbstractDiscretization{T,N},A<:AbstractArray{T,N}}

Produces a discretization of the given variable based on `T` and array type `A`.
"""
discretize(::Type{A}, ::D, ::Var) where {Q,T,N,D<:AbstractDiscretization{Q,N},A<:AbstractArray{T,N}} = error("missing discretize implementation for $D")
discretize(d::AbstractDiscretization{Q,N}, var::Var) where {Q,N} = discretize(Array{vartype(var),N}, d, var)
discretize(::Type{A}, grid::Grid, var::Var) where {A<:AbstractVector} = zero(similar(A{vartype(var)}, dimlength(var.dim, grid)))

# prognostic state vector constructor
function prognosticstate(::Type{A}, grid::Grid, layervars::NamedTuple, gridvars::Tuple) where {T,A<:AbstractArray{T}}
    # get lengths
    gridvar_sizes = map(v -> dimlength(vardims(v), grid), gridvars)
    layervar_sizes = map(vars -> map(v -> dimlength(vardims(v), grid), vars), layervars)
    Ng = length(gridvar_sizes) > 0 ? sum(gridvar_sizes) : 0
    Nl = sum(map(vars -> length(vars) > 0 ? sum(vars) : 0, layervar_sizes))
    # build axis indices;
    # non-grid prognostic variables get collected at the top of the vector, in the order provided
    i = 1
    layervar_ax = map(layervars, layervar_sizes) do vars, sizes
        j = 1
        coords = map(vars, sizes) do var, N
            coord = varname(var) => j:j+N-1
            j += N
            return coord
        end
        i += j
        return (; coords...)
    end
    # grid variables get interlaced throughout the rest of the vector; i.e. for variable i, its grid points are:
    # i:k:kn where k is the number of grid variables and n is the length of the grid.
    gridvar_ax = (;(varname(p) => st:length(gridvars):(Ng+Nl) for (p,st) in zip(gridvars, (Nl+1):(1+Nl+length(gridvars))))...)
    # select only non-empty layers
    layervar_ax = (;(name => layervar_ax[name] for name in keys(layervar_ax) if length(layervar_ax[name]) > 0)...)
    # allocate component array; assumes all prognostic variables have the same type (and they should!)
    u = zero(similar(A, Ng+Nl))
    u_ax = map(ax -> ViewAxis(first(ax)[1]:last(ax)[end], Axis(ax)), layervar_ax)
    return ComponentVector(u, (Axis(merge(u_ax, gridvar_ax)),))
end
