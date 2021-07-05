abstract type GridSpec end
struct Edges <: GridSpec end
struct Cells <: GridSpec end

export GridSpec, Edges, Cells

abstract type Geometry end
struct UnitVolume <: Geometry end

"""
    struct Grid{S,G,A,Q} <: AbstractVector{Q}

Represents the 1D spatial discretization on which time integration is performed. `S` is a `GridSpec`,
either `Edges` or `Cells` (always edges upon initial construction). The grid representation can be
converted (allocation free) between grid edges and cells via the `cells` and `edges` methods. `G`
represents the geometry/volume on which the vertical 1D discretization is applied. `A` is the underlying
array type, and `Q` is the numerical type (e.g. `Float64` or a `Unitful.Quantity`).
"""
struct Grid{S,G,A,Q} <: AbstractVector{Q}
    geometry::G
    values::NamedTuple{(:edges,:cells),NTuple{2,SubArray{Q,1,A,Tuple{StepRange{Int64,Int64}},true}}}
    deltas::NamedTuple{(:edges,:cells),NTuple{2,SubArray{Q,1,A,Tuple{StepRange{Int64,Int64}},true}}}
    function Grid(vals::TArray, geometry::G=UnitVolume()) where {G<:Geometry,Q<:Number,TArray<:AbstractArray{Q,1}}
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
        new{Edges,G,typeof(edges),Q}(geometry,
            (edges=view(edges,1:1:nedges),cells=view(cells,1:1:ncells)),
            (edges=view(Δedges,1:1:nedges-1),cells=view(Δcells,1:1:ncells-1)),
        )
    end
    function Grid(grid::Grid{Edges,G,A,Q}, interval::ClosedInterval{Int}) where {G,A,Q}
        start, stop = interval.left, interval.right
        edges = @view grid.values.edges[start:stop]
        cells = @view grid.values.cells[start:stop-1]
        Δedges = @view grid.deltas.edges[start:stop-1]
        Δcells = @view grid.deltas.cells[start:stop-2]
        new{Edges,G,A,Q}(grid.geometry, (edges=edges,cells=cells), (edges=Δedges,cells=Δcells))
    end
    function Grid(grid::Grid{Cells,G,A,Q}, interval::ClosedInterval{Int}) where {G,A,Q}
        start, stop = interval.left, interval.right
        edges = @view grid.values.edges[start:stop+1]
        cells = @view grid.values.cells[start:stop]
        Δedges = @view grid.deltas.edges[start:stop]
        Δcells = @view grid.deltas.cells[start:stop-1]
        new{Cells,G,A,Q}(grid.geometry, (edges=edges,cells=cells), (edges=Δedges,cells=Δcells))
    end
    Grid(::Type{Cells}, grid::Grid{Edges,G,A,Q}) where {G,A,Q} = new{Cells,G,A,Q}(grid.geometry,grid.values,grid.deltas)
    Grid(::Type{Edges}, grid::Grid{Cells,G,A,Q}) where {G,A,Q} = new{Edges,G,A,Q}(grid.geometry,grid.values,grid.deltas)
end

Base.show(io::IO, grid::Grid{S,G}) where {S,G} = print(io, "Grid{$S}($(grid[1])..$(grid[end])) of length $(length(grid)) with geometry $G")
Base.show(io::IO, ::MIME{Symbol("text/plain")}, grid::Grid) = show(io, grid)

function subgrid(grid::Grid{S,G,A,Q}, interval::Interval{L,R,Q}) where {S,G,A,Q,L,R}
    l,r = interval.left,interval.right
    vals = values(grid)
    # Determine indices which lie in the given interval
    l_ind = findfirst(x -> x ∈ interval, vals)
    r_ind = findlast(x -> x ∈ interval, vals)
    @assert !isnothing(l_ind) && !isnothing(r_ind) "No grid points in the given interval $interval"
    # Map back to full grid indices
    Grid(grid,l_ind..r_ind)
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
@propagate_inbounds Base.getindex(grid::Grid{S,G,A,Q}, interval::Interval{L,R,Q}) where {S,G,A,Q,L,R} = subgrid(grid,interval)
Base.setindex!(grid::Grid, args...) = error("setindex! is not allowed for Grid types")

regrid(x::AbstractVector, xgrid::Grid, newgrid::Grid, interp=Linear(), bc=Line()) =
    regrid!(similar(x,length(newgrid)), x, xgrid, newgrid, interp, bc)
function regrid!(out::AbstractVector, x::AbstractVector, xgrid::Grid, newgrid::Grid, interp=Linear(), bc=Line())
    let f = @> interpolate((xgrid,), x, Gridded(interp)) extrapolate(bc);
        out .= f.(newgrid)
    end
end

export Grid, cells, edges, indexmap, subgrid, Δ, regrid, regrid!

"""
    Profile(pairs...;names)

Constructs a Profile from the given pairs Q => (x1,...,xn) where x1...xn are the values defined at Q.
Column names for the resulting DimArray can be set via the names parameter which accepts an NTuple of symbols,
where N must match the number of parameters given (i.e. n).
"""
function Profile(pairs::Pair{Q,NTuple{N,T}}...;names::Union{Nothing,NTuple{N,Symbol}}=nothing) where {T,N,Q<:DistQuantity}
    depths, vals = zip(pairs...)
    params = hcat(collect.(vals)...)'
    names = isnothing(names) ? [Symbol(:x,:($i)) for i in 1:N] : collect(names)
    DimArray(params, (depth=collect(depths), var=names))
end

"""
    interpolateprofile(profile::Profile, state; interp=Linear())

Interpolates the given profile to the corresponding variable grids. Assumes state to be indexable via the corresponding
variable symbol and that the parameter names in state and profile match.
"""
function interpolateprofile!(profile::DimArray, state; interp=Linear())
    let (depths,names) = dims(profile),
        z = ustrip.(depths);
        for p in names
            f = @> interpolate((z,), profile[:,p], Gridded(interp)) extrapolate(Flat())
            state[p] .= f.(state.grids[p]) .|> dustrip   # assume length(grid) == length(state.p)
        end
    end
end

export Profile, interpolateprofile!
