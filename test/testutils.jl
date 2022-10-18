using CryoGrid.Numerics
using CryoGrid.Utils

import Unitful

# workaround for @test not allowing error messages
printerror(str) = begin (@error str); false; end

function allfinite(x::AbstractArray)
	chk = isfinite.(x)
	all(chk) || printerror("Infinite values at $(findall(.!chk)) in $x")
end

function allequal(x::AbstractArray, y::AbstractArray; atol=1.0e-3)
	chk = .≈(x,y,atol=atol)
	all(chk) || printerror("Values do not match at $(findall(.!chk)) in $x and $y")
end

function build_test_state(grid::Grid, layer::Layer)
	vargrid(::OnGrid{Cells}, grid::Grid) = cells(grid)
	vargrid(::OnGrid{Edges}, grid::Grid) = edges(grid)
	named_layer = Named(:layer, layer)
	vars = Strat._collectvars(named_layer)
	return (
		grid = grid,
		grids = (; map(v -> varname(v) => vargrid(vardims(v), grid), filter(isongrid, vars))...),
		map(v -> varname(v) => zeros(dimlength(vardims(v), grid))*varunits(v), vars)...
	)
end
