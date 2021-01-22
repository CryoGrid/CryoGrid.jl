using Base: @propagate_inbounds

abstract type GridSpec end
struct Edges <: GridSpec end
struct Cells <: GridSpec end

const CAxis = ComponentArrays.Axis

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
        deltas = similar(fullgrid, N-2)
        dedges = @view deltas[1:2:end]
        dedges .= edges[2:end] .- edges[1:end-1]
        dcells = @view deltas[2:2:end-1]
        dcells .= cells[2:end] .- cells[1:end-1]
        fullgrid = fullgrid |> SVector{N,Q}
        deltas = deltas |> SVector{N-2,Q}
        values = ComponentArray(fullgrid,(CAxis{(edges=1:2:N,cells=2:2:N-1)}(),))
        deltas = ComponentArray(deltas,(CAxis{(edges=1:2:N-2,cells=2:2:N-3)}(),))
        new{Edges,Q,typeof(values),typeof(deltas)}(values,deltas)
    end
    Grid(::Type{Cells}, grid::Grid{Edges,Q,G,D}) where {Q,G,D} =
        new{Cells,Q,G,D}(grid.values,grid.deltas)
    Grid(::Type{Edges}, grid::Grid{Cells,Q,G,D}) where {Q,G,D} =
        new{Edges,Q,G,D}(grid.values,grid.deltas)
end

function findinterval(values::A, interval::ClosedInterval{Q}) where {Q,A<:AbstractArray{Q}}
    l,r = interval.left,interval.right
    # Determine shortest interval that encloses l and r in *full* grid
    l_ind = argmin(abs.(values .- l))+1
    r_ind = argmin(abs.(values .- r))
    values[l_ind:r_ind]
end
# TODO: there must be a better way to do this than repeating defs for both edges and cells;
# Maybe a macro?
Δ(grid::Grid{Edges}) = grid.deltas.edges
Δ(grid::Grid{Cells}) = grid.deltas.cells
cells(grid::Grid{Edges}) = Grid(Cells, grid)
edges(grid::Grid{Cells}) = Grid(Edges, grid)
subgrid(grid::Grid{Edges,Q}, interval::ClosedInterval{Q}) where Q = Grid(findinterval(grid.values.edges, interval))
subgrid(grid::Grid{Cells,Q}, interval::ClosedInterval{Q}) where Q = Grid(findinterval(grid.values.cells, interval))
Base.similar(grid::Grid{Edges}) = Grid(grid.values.edges)
Base.similar(grid::Grid{Cells}) = Grid(grid.values.edges) |> cells
Base.size(grid::Grid{Edges}) = size(grid.values.edges)
Base.size(grid::Grid{Cells}) = size(grid.values.cells)
Base.length(grid::Grid{Edges}) = length(grid.values.edges)
Base.length(grid::Grid{Cells}) = length(grid.values.cells)
@propagate_inbounds Base.getindex(grid::Grid{Edges}, i::Int) = grid.values.edges[i]
@propagate_inbounds Base.getindex(grid::Grid{Cells}, i::Int) = grid.values.cells[i]
@propagate_inbounds Base.getindex(grid::Grid{S,Q}, interval::ClosedInterval{Q}) where {S,Q} = subgrid(grid,interval)
AxisArrays.axistrait(::Type{<:Grid}) = AxisArrays.Dimensional
AxisArrays.axisindexes(::Type{AxisArrays.Dimensional}, ax::Grid{Edges}, idx) =
    AxisArrays.axisindexes(AxisArrays.Dimensional,ax.values.edges,idx)
AxisArrays.axisindexes(::Type{AxisArrays.Dimensional}, ax::Grid{Cells}, idx) =
    AxisArrays.axisindexes(AxisArrays.Dimensional,ax.values.cells,idx)

export GridSpec, Edges, Cells
export Grid, cells, edges, Δ
