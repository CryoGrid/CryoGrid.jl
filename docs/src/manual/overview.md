# [Overview](@id overview)
## Setting up a model

```@meta
DocTestSetup = quote
    using CryoGrid
end
```

At the highest level, a model in `CryoGrid.jl` is defined by one or more [`Tile`](@ref)s each consisting of a [`Grid`](@ref) and a [`Stratigraphy`](@ref), constructed top-down from individual [`Layer`](@ref)s, each of which has one or more [`Process`](@ref)es. Each layer in the `Stratigraphy` is assigned a depth, which then aligns it with the `Grid`. All models must consist of at least three layers/nodes: `Top` and `Bottom` layers with corresponding boundary conditions, as well as one or more [`SubSurface`](@ref) layers. Here we define a simple three-layer model (or one-layer, exlcuding the boundaries) with a single sub-surface process, i.e. [`HeatBalance`](@ref) (heat conduction):

```julia
# ... load forcings, set up profiles, etc.
# see examples/heat_vgfc_seb_saoylov_custom.jl for more details
forcings = (;Tair,pr,q,wind,Lin,Sin,z)
strat = Stratigraphy(
    -2.0u"m" => Top(SurfaceEnergyBalance(forcings)),
    0.0u"m" => Ground(soilprofile, HeatBalance(:H; freezecurve=DallAmico())),
    1000.0u"m" => Bottom(GeothermalHeatFlux(0.053u"J/s/m^2"))
);
grid = CryoGrid.DefaultGrid_5cm
# define initial conditions for temperature using a given profile;
# The default initializer linearly interpolates between profile points.
initT = initializer(:T, tempprofile)
tile = Tile(strat, grid, initT);
```

This model can then be used to construct a `CryoGridProblem`:

```julia
tspan = (DateTime(2010,10,30),DateTime(2011,10,30))
p = parameters(tile)
u0 = initialcondition!(tile, tspan, p, initT)
prob = CryoGridProblem(tile, u0, tspan, p, saveat=24*3600.0, savevars=(:T,)) # produces an ODEProblem with problem type CryoGridODEProblem
```

It can then be solved/integrated using the `solve` function (from `DiffEqBase` and `OrdinaryDiffEq`):

```julia
# solve and construct CryoGridOutput from solution
sol = @time solve(prob, saveat=24*3600.0, progress=true);
out = CryoGridOutput(sol)
```

The resulting `CryoGridOutput` type provides `DimArray`s containing the model outputs over time and space:

```raw
julia> out.T
278×366 DimArray{Float64,2} with dimensions: 
  Z: Quantity{Float64, 𝐋, Unitful.FreeUnits{(m,), 𝐋, nothing}}[0.01 m, 0.03 m, …, 850.0 m, 950.0 m] Sampled: Ordered Irregular Points,
  Ti (Time): DateTime[2010-10-30T00:00:00, …, 2011-10-30T00:00:00] Sampled: Ordered Irregular Points
```

## Defining model behavior

Notice that, in the example above, it is types such as `Ground`, `HeatBalance`, `DallAmico`, etc. that specify which components the model should use. These components are defined by adding method dispatches to the [CryoGrid interface](@ref toplevel) methods. State variables are declared via the [`variables`](@ref) method, e.g:

```julia
variables(soil::Soil, heat::HeatBalance{<:EnthalpyBased}) = (
    Prognostic(:H, OnGrid(Cells), u"J/m^3"),
    Diagnostic(:T, OnGrid(Cells), u"°C"),
    Diagnostic(:C, OnGrid(Cells), u"J//K*/m^3"),
    Diagnostic(:∂H∂T, OnGrid(Cells), u"J/K/m^3"),
    Diagnostic(:k, OnGrid(Edges), u"W/m/K"),
    Diagnostic(:kc, OnGrid(Cells), u"W//m/K"),
)
```

When the `HeatBalance` process is assigned to a `Soil` layer, `Tile` will invoke this method and create state variables corresponding to each [`Var`](@ref). [`Prognostic`](@ref) variables are assigned derivatives (in this case, `dH`, since `H` is the prognostic state variable) and integrated over time. `Diagnostic` variables provide in-place caches for derived/intermediary state variables.

Each variable definition consists of a name (a Julia `Symbol`), a type, and a shape. For variables discretized on the grid, the shape is specified by `OnGrid`, which will generate an array of the appropriate size when the model is compiled. The arguments `Cells` and `Edges` specify whether the variable should be defined on the grid cells or edges respecitvely.

The real work finally happens in [`computediagnostic!`](@ref) and [`computeprognostic!`](@ref), the latter of which should be used to compute the time derivatives (here `dH`). [`interact!`](@ref) defines the behavior at the boundaries and should be used to compute the derivatives (and any other necessary values) at the interface between layers.

We can take as an example the implementation of `computeprognostic!` for enthalpy-based heat conduction (note that `jH` is a diagnostic variable representing the energy flux over each cell edge):

```julia
function CryoGrid.computeprognostic!(::SubSurface, ::HeatBalance{<:EnthalpyBased}, state)
    Δk = Δ(state.grid) # cell sizes
    ΔT = Δ(cells(state.grid)) # midpoint distances
    # compute internal fluxes and non-linear diffusion assuming boundary fluxes have been set
    nonlineardiffusion!(state.dH, state.jH, state.T, ΔT, state.k, Δk)
    return nothing
end
```

!!! warning

    Prognostic state variables like `H` in the example above **should not be directly modified** in the model code. They should only be modified by the calling solver/integrator. This is especially important when using higher order or implicit integrators as unexpected changes to prognostic state may destroy the accuracy of their internal interpolant. For modeling discontinuities, use [`Events`](@ref) instead.

Note that `state` is (typically) of type [`LayerState`](@ref) with properties corresponding to the state variables declared by the `variables` function for `Soil` and `HeatBalance`. Additionally, output arrays for the time derivatives are provided (here `dH`), as well as the current timestep, layer boundary depths, and variable grids (accessible via `state.t`, and `state.grid` respectively). Note that `state` will also contain other variables declared on this `Soil` layer by other `SubSurfaceProcess`es, allowing for implicit coupling between processes where appropriate.
