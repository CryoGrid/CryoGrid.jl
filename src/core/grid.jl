using Base: @propagate_inbounds

abstract type GridSpec end
struct Edges <: GridSpec end
struct Cells <: GridSpec end

export GridSpec, Edges, Cells

struct Grid{S,Q,G,D} <: AbstractArray{Q,1}
    values::G
    deltas::D
    Grid(vals::TArray) where {Q,TArray<:AbstractArray{Q,1}} = begin
        # initialize full grid (edges with cells/midpoints)
        N = 2*length(vals)-1
        fullgrid = MVector{N,Q}(undef)
        edges = @view fullgrid[1:2:end]
        edges .= vals
        cells = @view fullgrid[2:2:end-1]
        cells .= (vals[1:end-1] + vals[2:end]) / (2.0*one(Q))
        # initialize delta grid; note that the distinction between edges and cells here is kind of
        # redundant since the distance between midpoints should always be the same as distance between boundaries...
        deltas = similar(fullgrid, N-2)
        dedges = @view deltas[1:2:end]
        dedges .= edges[2:end] .- edges[1:end-1]
        dcells = @view deltas[2:2:end-1]
        dcells .= cells[2:end] .- cells[1:end-1]
        # convert to immutable SVectors
        fullgrid = fullgrid |> SVector{N,Q}
        deltas = deltas |> SVector{N-2,Q}
        # manually specify interleaved axes
        values = ComponentArray(fullgrid,(CAxis{(edges=1:2:N,cells=2:2:N-1)}(),))
        deltas = ComponentArray(deltas,(CAxis{(edges=1:2:N-2,cells=2:2:N-3)}(),))
        new{Edges,Q,typeof(values),typeof(deltas)}(values,deltas)
    end
    Grid(grid::Grid{S,Q,G,D}, inds::AbstractRange) where {S<:GridSpec,Q,G,D} = begin
        start, stop = inds.start, inds.stop
        # create new ComponentArrays with adjusted axes; this should be allocation-free
        values = ComponentArray(grid.values,(CAxis{(edges=start:2:stop,cells=start+1:2:stop-1)}(),))
        deltas = ComponentArray(grid.deltas,(CAxis{(edges=start:2:stop-2,cells=start+1:2:stop-3)}(),))
        new{S,Q,typeof(values),typeof(deltas)}(values,deltas)
    end
    Grid(::Type{Cells}, grid::Grid{Edges,Q,G,D}) where {Q,G,D} =
        new{Cells,Q,G,D}(grid.values,grid.deltas)
    Grid(::Type{Edges}, grid::Grid{Cells,Q,G,D}) where {Q,G,D} =
        new{Edges,Q,G,D}(grid.values,grid.deltas)
end

function subgrid(grid::Grid{S,Q}, interval::ClosedInterval{Q}) where {S,Q}
    l,r = interval.left,interval.right
    vals = values(grid)
    # Determine shortest interval that encloses l and r in *full* grid
    l_ind = argmin(abs.(vals .- l))+1
    r_ind = argmin(abs.(vals .- r))
    # Map back to full grid indices
    idxmap = indexmap(grid)
    Grid(grid,idxmap[l_ind]:idxmap[r_ind])
end
Δ(grid::Grid{Edges}) = grid.deltas.edges
Δ(grid::Grid{Cells}) = grid.deltas.cells
cells(grid::Grid{Edges}) = Grid(Cells, grid)
edges(grid::Grid{Cells}) = Grid(Edges, grid)
@inline indexmap(grid::Grid{Edges}) = ComponentArrays.indexmap(getaxes(grid.values)[1]).edges
@inline indexmap(grid::Grid{Cells}) = ComponentArrays.indexmap(getaxes(grid.values)[1]).cells
@inline Base.values(grid::Grid{Edges}) = grid.values.edges
@inline Base.values(grid::Grid{Cells}) = grid.values.cells
@inline Base.similar(grid::Grid{Edges}) = Grid(grid.values.edges)
@inline Base.similar(grid::Grid{Cells}) = Grid(grid.values.edges) |> cells
@inline Base.size(grid::Grid) = size(values(grid))
@inline Base.length(grid::Grid) = length(values(grid))
@propagate_inbounds Base.getindex(grid::Grid, i::Int) = values(grid)[i]
@propagate_inbounds Base.getindex(grid::Grid{S,Q}, interval::ClosedInterval{Q}) where {S,Q} = subgrid(grid,interval)
AxisArrays.axistrait(::Type{<:Grid}) = AxisArrays.Dimensional
AxisArrays.axisindexes(::Type{AxisArrays.Dimensional}, ax::Grid, idx) =
    AxisArrays.axisindexes(AxisArrays.Dimensional,values(ax),idx)

export Grid, cells, edges, indexmap, subgrid, Δ
