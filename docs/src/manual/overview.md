# Overview
## Setting up a model

At the highest level, a model in `CryoGrid.jl` is defined by a [`Grid`](@ref) and a [`Stratigraphy`](@ref), constructed top-down from individual [`StratNode`](@ref)s, each of which consists of a [`Layer`](@ref) and one or more [`Process`](@ref)es. Each node in the `Stratigraphy` is assigned a depth, which then aligns it with the `Grid`. All models must consist of at least three layers/nodes: `Top` and `Bottom` layers with corresponding boundary conditions, as well as one or more [`SubSurface`](@ref) layers. Here we define a simple three-layer model (or one-layer, exlcuding the boundaries) with a single sub-surface process, i.e. [`Heat`](@ref) (heat conduction):

```julia
# ... load forcings, set up profiles, etc.
# see examples/heat_vgfc_seb_saoylov_custom.jl for more details
strat = Stratigraphy(
    -2.0u"m" => Top(SurfaceEnergyBalance(Tair,pr,q,wind,Lin,Sin,z)),
    0.0u"m" => Ground(:soil, Soil(soilprofile), Heat{:H}(tempprofile, freezecurve=SFCC(VanGenuchten()))),
    1000.0u"m" => Bottom(GeothermalHeatFlux(0.053u"J/s/m^2"))
);
grid = CryoGrid.Models.DefaultGrid_5cm
model = CryoGridSetup(strat,grid);
```

This model can then be used to construct an `ODEProblem` (from `DiffEqBase.jl`) via the `CryoGridProblem` constructor:

```julia
tspan = (DateTime(2010,10,30),DateTime(2011,10,30))
prob = CryoGridProblem(model,tspan) # produces an ODEProblem with problem type CryoGridODEProblem
```

It can then be solved simply using the `solve` function (also from `DiffEqBase` and `OrdinaryDiffEq`):

```julia
# solve with forward Euler (fixed 5 minute time steps) and construct CryoGridOutput from solution
out = @time solve(prob, Euler(), dt=5*60.0, saveat=24*3600.0, progress=true) |> CryoGridOutput;
```

The result is a `CryoGridOutput` type which provides access to `DimArrays` containing the model outputs over time and space:

```
julia> out.soil.T
278×366 DimArray{Float64,2} with dimensions: 
  Z: Quantity{Float64, 𝐋, Unitful.FreeUnits{(m,), 𝐋, nothing}}[0.01 m, 0.03 m, …, 850.0 m, 950.0 m] Sampled: Ordered Irregular Points,
  Ti (Time): DateTime[2010-10-30T00:00:00, …, 2011-10-30T00:00:00] Sampled: Ordered Irregular Points
```

## Defining model behavior

Notice that, in the example above, it is types such as `Soil`, `Heat`, `SFCC`, etc. that specify which components the model should use. These components are defined by adding method dispatches to the functions defined in the [`Interface`](@ref) module. State variables are declared via the [`variables`](@ref) method, e.g:

```julia
""" Variable definitions for heat conduction (enthalpy) on soil layer. """
variables(soil::Soil, heat::Heat{:H}) = (
    Prognostic(:H, Float"J/m^3", OnGrid(Cells)),
    Diagnostic(:T, Float"K", OnGrid(Cells)),
    Diagnostic(:C, Float"J//K*/m^3", OnGrid(Cells)),
    Diagnostic(:Ceff, Float"J/K/m^3", OnGrid(Cells)),
    Diagnostic(:k, Float"W/m/K", OnGrid(Edges)),
    Diagnostic(:kc, Float"W//m/K", OnGrid(Cells)),
    # this last line just appends any state variables or parameters
    # defined by the freeze curve to the tuple.
    variables(soil, heat, freezecurve(heat))...,
)
```

When the `Heat` process is assigned to a `Soil` layer, `CryoGridSetup` will invoke this method and create state variables corresponding to each [`Var`](@ref). [`Prognostic`](@ref) variables are assigned derivatives (in this case, `dH`, since `H` is the prognostic state variable) and integrated over time. `Diagnostic` variables provide in-place caches for intermediary variables/computations and are automatically tracked by the modeling engine (i.e. their saved values will appear in `CryoGridOutput`).

Each variable definition consists of a name (a Julia `Symbol`), a type, and a shape. For variables discretized on the grid, the shape is specified by `OnGrid`, which will generate an array of the appropriate size when the model is compiled. The arguments `Cells` and `Edges` specify whether the variable should be defined on the grid cells or edges respecitvely.

The real work finally happens in [`diagnosticstep!`](@ref) and [`prognosticstep!`](@ref), the latter of which should be used to compute the time derivatives (here `dH`). [`interact!`](@ref) defines the behavior at the boundaries and should be used to compute the derivatives (and any other necessary values) at the interface between layers.

We can take as an example the implementation of `prognosticstep!` for enthalpy-based heat conduction:

```julia
""" Prognostic step for heat conduction (enthalpy) on soil layer. """
function prognosticstep!(::Soil, ::Heat{:H}, state)
    Δk = Δ(state.grids.k) # cell sizes
    ΔT = Δ(state.grids.T)
    # Diffusion on non-boundary cells
    heatconduction!(state.dH,state.T,ΔT,state.k,Δk)
end
```

!!! warning

    Prognostic state variables like `H` in the example above **should not be directly modified** in user code. This is especially crucial when using higher order or implicit integrators as unexpected changes to the underlying state may destroy the accuracy of their internal interpolators. For modleing discontinuities, use [`Callbacks`](@ref) instead.

Note that `state` is a `NamedTuple` with fields corresponding to the variables declared by the `variables` function for `Soil` and `Heat`. Additionally, output arrays for the time derivatives are provided (here `dH`), as well as the current timestep, layer upper boundary depth, parameters, and variable grids (accessible via `state.t`, `state.z`, `state.params`, and `state.grids` respectively). Note that `state` will also contain other variables declared on this `Soil` layer by other `SubSurfaceProcess`es, allowing for implicit coupling between processes where appropriate.
