# Overview
## Setting up a model

At the highest level, a model in `CryoGrid.jl` is defined by a [`Grid`](@ref) and a [`Stratigraphy`](@ref), constructed top-down from individual [`StratNode`](@ref)s, each of which consists of a [`Layer`](@ref) and one or more ['Process`](@ref)es. Each node in the `Stratigraphy` is assigned a depth, which then aligns it with the `Grid`. All models must consist of at least three layers/nodes: `Top` and `Bottom` layers with corresponding boundary conditions, as well as one or more [`SubSurface`](@ref) layers. Here we define a simple three-layer model (or one-layer, exlcuding the boundaries) with a single sub-surface process, i.e. [`Heat`](@ref) (heat conduction):

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

In the example above, we are using types such as `Soil`, `Heat`, `SFCC`, etc. to specify which components the model should use. These components are defined by adding method dispatches to the functions defined in the [`Interface`](@ref) module. State variables are declared via the [`variables`](@ref) method, e.g:

```julia
""" Variable definitions for heat conduction (enthalpy) on soil layer. """
variables(soil::Soil, heat::Heat{:H}) = (
    Prognostic(:H, Float"J/m^3", OnGrid(Cells)),
    Diagnostic(:T, Float"K", OnGrid(Cells)),
    Diagnostic(:C, Float"J//K*/m^3", OnGrid(Cells)),
    Diagnostic(:Ceff, Float"J/K/m^3", OnGrid(Cells)),
    Diagnostic(:k, Float"W/m/K", OnGrid(Edges)),
    Diagnostic(:kc, Float"W//m/K", OnGrid(Cells)),
    variables(soil, heat, freezecurve(heat))...,
)
```

When the `Heat` process is assigned to a `Soil` layer, `CryoGridSetup` will invoke this method and create state variables corresponding to each [`Var`](@ref). [`Prognostic`](@ref) variables are assigned derivatives (in this case, `dH`, since `H` is the prognostic state variable) and integrated over time. `Diagnostic` variables provide in-place caches for intermediary variables/computations and are automatically tracked by the modeling engine (i.e. their saved values will appear in `CryoGridOutput`).

Computations are performed in [`diagnosticstep!`](@ref) and [`prognosticstep!`](@ref), the latter of which should be used to compute the time derivatives (here `dH`). [`interact!`](@ref) defines the behavior at the boundaries and should be used to compute the derivatives (and any other necessary values) at the edges of each layer.
