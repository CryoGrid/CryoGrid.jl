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

function build_test_state(grid::Grid, layer::Layer, name::Symbol=:layer; t=0.0, with_units=true)
	vargrid(::OnGrid{Cells}, grid::Grid) = cells(grid)
	vargrid(::OnGrid{Edges}, grid::Grid) = edges(grid)
	maybeunits(var::Var) = with_units ? varunits(var) : Unitful.NoUnits
	named_layer = Named(name, layer)
	vars = CryoGrid.variables(named_layer)
	return (
		t = t,
		grid = grid,
		grids = (; map(v -> varname(v) => vargrid(vardims(v), grid), filter(isongrid, vars))...),
		map(v -> varname(v) => zeros(dimlength(vardims(v), length(grid)))*maybeunits(v), vars)...
	)
end
