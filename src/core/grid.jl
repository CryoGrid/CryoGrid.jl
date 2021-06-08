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
        fullgrid = similar(vals,N)
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
        # fullgrid = fullgrid |> SVector{N,Q}
        # deltas = deltas |> SVector{N-2,Q}
        # manually specify interleaved axes
        values = ComponentArray(fullgrid,(Axis{(edges=1:2:N,cells=2:2:N-1)}(),))
        deltas = ComponentArray(deltas,(Axis{(edges=1:2:N-2,cells=2:2:N-3)}(),))
        new{Edges,Q,typeof(values),typeof(deltas)}(values,deltas)
    end
    Grid(grid::Grid{S,Q,G,D}, interval::ClosedInterval{Int}) where {S<:GridSpec,Q,G,D} = begin
        start, stop = interval.left, interval.right
        # create new ComponentArrays with adjusted axes; this should be allocation-free
        let values = getdata(grid.values),
            deltas = getdata(grid.deltas),
            valaxis = Axis((edges=start:2:stop,cells=start+1:2:stop-1)),
            delaxis = Axis((edges=start:2:stop-2,cells=start+1:2:stop-3)),
            newvalues = ComponentArray(values,(valaxis,)),
            newdeltas = ComponentArray(deltas,(delaxis,));
            new{S,Q,typeof(newvalues),typeof(newdeltas)}(newvalues,newdeltas)
        end
    end
    Grid(::Type{Cells}, grid::Grid{Edges,Q,G,D}) where {Q,G,D} =
        new{Cells,Q,G,D}(grid.values,grid.deltas)
    Grid(::Type{Edges}, grid::Grid{Cells,Q,G,D}) where {Q,G,D} =
        new{Edges,Q,G,D}(grid.values,grid.deltas)
end

Base.show(io::IO, grid::Grid{S}) where S = print(io, "Grid{$S}($(grid[1])..$(grid[end])) of length $(length(grid))")
Base.show(io::IO, mime::MIME{Symbol("text/plain")}, grid::Grid{S}) where S = show(io, grid)
Base.show(io::IO, ::Type{<:Grid{S,T}}) where {S,T} = print(io, "Grid{$S,$T,...}")

function subgrid(grid::Grid{S,Q}, interval::Interval{L,R,Q}) where {S,L,R,Q}
    l,r = interval.left,interval.right
    vals = values(grid)
    # Determine indices which lie in the given interval
    l_ind = findfirst(x -> x ∈ interval, vals)
    r_ind = findlast(x -> x ∈ interval, vals)
    @assert !isnothing(l_ind) && !isnothing(r_ind) "No grid points in the given interval $interval"
    # Map back to full grid indices
    idxmap = indexmap(grid)
    Grid(grid,idxmap[l_ind]..idxmap[r_ind])
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
@propagate_inbounds Base.getindex(grid::Grid{S,Q}, interval::Interval{L,R,Q}) where {S,L,R,Q} = subgrid(grid,interval)

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
    DimArray(params, (Z(collect(depths)), Y(names)))
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
