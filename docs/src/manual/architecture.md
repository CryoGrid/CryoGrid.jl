# [Architecture](@id arch)

This page provides a general overview of the code organization and architecture of CryoGrid.jl.

## Modular design

[Modular programming](https://en.wikipedia.org/wiki/Modular_programming) in software design revolves around the separation of computer programs into "modules" that can be independently constructed, tested, and coupled together into a larger system. The benefits of modular programming are well documented and have been common practice in software engineering for decades [1]. In the context of physical modeling, modular programming has the potential to facilitate rapid prototyping and comparison of different model configurations, parameterizations, and process interactions [2].

The term "modular programming" is fairly abstract and encompasses a wide range of patterns and techniques centered around the central aim of building robust and reusable software components. [Namespaces](https://en.wikipedia.org/wiki/Namespace) are a commonly employed tool for organizing code into standalone modules or packages that share functionality and naming patterns. Namespaces help avoid name collisions by localizing variable and function names to their enclosing namespace, therefore resolving possible ambiguities.

In Julia, namespaces are declared via [modules](https://docs.julialang.org/en/v1/manual/modules/). Modules are more-or-less self-contained namespaces which can be used to organize and isolate code. Modules can "export" methods or variables which are intended for external use via `export` statements. Other modules can then import these methods or variables into their namespace via `using` and `import` statements; e.g. `using Dates` brings all `export`ed names from the `Dates` module into the current namespace. Note that the top-level module (i.e. in a script or in the REPL) is always called `Main`.

The `CryoGrid` module provided by CryoGrid.jl is organized into a series of submodules:

|Name|Description|Depends on|
|----|-----------|----------|
|`Utils`|Miscellaneous utility methods and types.||
|`Numerics`|Utilities for math, array caches, and spatial discretization.|`Utils`|
|`InputOutput`|Methods and types related to reading and writing input and output data.|`Utils`,`Numerics`|
|`Diagnostics`|Tools for model diagnostics.|`Utils`,`Numerics`|
|`Hydrology`|Methods and types for computing water related quanities.|`Utils`,`Numerics`|
|`Heat`|Methods and types for computing heat and energy related quanities.|`Utils`,`Numerics`,`Hydrology`|
|`Soils`|Defines `Soil` layers and provides dispatches for soil-specific physical processes.|`Utils`,`Numerics`,`Hydrology`,`Heat`|
|`Snow`|Defines `Snowpack` layer and provides dispatches for snow processes.|`Utils`,`Numerics`,`Hydrology`,`Heat`|
|`Salt`|Provides types and dispatches for coupled heat/salt diffusion in saline soils.|`Utils`,`Numerics`,`Hydrology`,`Heat`,`Soils`|
|`Surface`|Defines boundary processes for the surface such as the surface energy and water balance equations.|`Utils`,`Numerics`,`Hydrology`,`Heat`,`Soils`,`Snow`|
|`Tiles`|Defines the `Tile` and `Stratigraphy` types for constructing 1D land models.|`Utils`,`Numerics`,`InputOutput`|
|`DiffEq`|Provides dispatches and utilities for integrating with solvers from the SciML `OrdinaryDiffEq` package.|`Utils`,`Numerics`,`InputOutput`|
|`LiteImplicit`|Provides an implementation of the `CryoGridLite` solver scheme from Langer et al. 2023.|`Utils`,`Numerics`|
|`Presets`|Provides pre-defined stratigraphies, forcings, and layer configurations to facilitate rapid prototyping.|`Utils`,`Numerics`,`InputOutput`,`Heat`,`Hydrology`,`Soils`|

Note that all submodules depend on the top-level `CryoGrid` module which declares all of the "core" types and [method interfaces](https://docs.julialang.org/en/v1/manual/interfaces/) for defining model behavior. Each submodule may additionally define its own method interfaces related to its own specific process(es) or layer(s).

The `@reexport` macro from the `Reexport` package is used extensively to propagate exported methods and types to the top-level `CryoGrid` namespace. This is intended to alleviate the user of the burden to keep track of which types/methods are exported by which submodules. In most cases, it is sufficient to simply `import` or `using` the `CryoGrid` module in order to bring all CryoGrid-related methods and types into scope.

## Model structure

In the context of CryoGrid, a "model" typically refers to one or more [`Tile`](@ref)s [2] which may or may not be laterally coupled together. A `Tile` typically corresponds to a rectangular volume discretized along the vertical z-axis, i.e. corresponding physically to depth/elevation. CryoGrid.jl implements this concept by defining a single `Tile` as a composition of the following:

- A [`Stratigraphy`](@ref) with three or more [`Layer`](@ref)s, including a [`Top`](@ref) layer and a [`Bottom`](@ref) layer.
- A [`StateVars`](@ref) cache which stores all non-prognostic state and grid data.
- Zero or more [`VarInitializer`](@ref)s that define the intial condition of the prognostic state.
- Zero or more layer [`Event`](@ref)s that may or may not be invoked when their trigger conditions are met.

The `Stratigraphy` is simply a `Tuple` of layers in ascending order of depth (i.e. top to bottom) paired with (initial) upper boundary depths. The thickness of each stratigraphy layer is therefore determined by the distance between the upper boundary of the layer and the upper boundary of the following layer. Depending on the configuration of the layer, this thickness may be either static or dynamic over time. In the latter case, the layer thickness `Δz` is automatically included as a prognostic state variable.

Each `SubSurface` layer in the stratigraphy will typically consist of one or more `Process`es as fields on the layer `struct` which should then be explicitly declared via a dispatch of the [`processes`](@ref) method. The [`variables`](@ref) and [`events`](@ref) methods similarly declare state variables and events respectively that should be defined for any given configuration of the layer.

The `Tile` constructor collects all of the relevant state variables declared by `variables` and discretizes them according to the given [`DiscretizationStrategy`](@ref). The resulting state vectors are initialized in the forward-diff compatible [`StateVars`](@ref) cache. On each invocation of `Tile`, the current [`TileState`](@ref) is constructed from the current prognostic state variable `u`, parameter vector `p`, and time step `t`. The `TileState` consists of named [`LayerState`](@ref)s which provide layer-local `view`s of each state variable array, i.e. only grid cells within the layer boundaries are included.

## Control flow

The `CryoGrid` module defines three primary methods that can be used to implement the behavior of each `Layer`/`Process` in any given model configuration. When updating the model at a particular timestep, these methods are typically invoked in the following order:

1. [`computediagnostic!`](@ref) updates all (non-flux) state variables and/or derived quantities based on the current (prognostic) state.
2. [`interact!`](@ref) defines interactions between adjacent layers in the stratigraphy, including fluxes over the layer boundary.
3. [`computeprognostic!`](@ref) computes all internal fluxes (and the divergence thereof) within each layer, after boundary fluxes are taken into account by `interact!`.

Layer and/or process specific implementations of each of these methods can generally assume that the previous methods have already been invoked by the caller (it is the responsibility of the calling code to ensure that this is the case). This is, for example, the order in which these methods will be invoked by `tile(du, u, p t)`.

Note that, due to the nature of multiple dispatch, the execution path (i.e. with respect to the actual source code) of any given model configuration will typically be quite nonlinear and may span multiple source files depending on where the matching method dispatches are defined. Users may find the `which` provided by Julia (and correspondingly the `@which` macro from `InteractiveUtils`) useful in figuring out which code is being executed. For example:

```julia
using CryoGrid
using CryoGrid.Diagnostics

soil = Ground()
grid = CryoGrid.Presets.DefaultGrid_5cm
state = Diagnostics.build_dummy_state(grid, soil)

@which CryoGrid.computediagnostic!(soil, state)
```
Output:
```
computediagnostic!(layer::Layer, state)
     @ CryoGrid ~/workspace/sparc-local/repos/CryoGrid/CryoGrid.jl/src/methods.jl:55
```

## State variables

In order to facilitate modularity and ease-of-use, CryoGrid.jl provides an automated system for initializing and configuring state variables for any given model configuration. Note that there is an important distinction between two types of model state: **prognostic** and **diagnostic**.

`Prognostic`(@ref) state variables fully define the state of the system at any given time `t`. They form what is typically called the "phase space" or "state space" in the mathematics and engineering literature. In order to be compatible with standard ODE solvers (e.g. like those in `OrdinaryDiffEq`), CryoGrid.jl automatically assembles prognostic state variables into a single array `u` (and its corresponding time derivative `du`) which is returned when initializing a `Tile` with the `initialcondition!` method. Note again that this array should always fully define the state of the system.

`Diagnostic`(@ref) state variables act as caches for intermediate and derived quantities defined by the model. They also may, in some cases, provide a means of coupling between different processes (e.g. the heat and water flux variables `jH` and `jw` might be updated by more than one `Process`). For any model configuration, all diagnostic variables should be fully updated (and thus consistent) with the given prognostic state after invoking `computediagnostic!`, `interact!`, and `computeprognostic!`.

When a `Tile` is constructed, all variables defined by each layer in the `Stratigraphy` are collected and then intiailized in [`StateVars`](@ref) according to the given `DiscretizationStrategy`.

## References

[1] Bass L, Clements P, Kazman R. Software architecture in practice. Addison-Wesley Professional; 2003.

[2] Westermann S, Ingeman-Nielsen T, Scheer J, Aalstad K, Aga J, Chaudhary N, Etzelmüller B, Filhol S, Kääb A, Renette C, Schmidt LS. The CryoGrid community model (version 1.0)–a multi-physics toolbox for climate-driven simulations in the terrestrial cryosphere. Geoscientific Model Development. 2023 May 15;16(9):2607-47.
