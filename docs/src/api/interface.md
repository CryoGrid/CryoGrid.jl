# Overview

At the highest level, a model in `CryoGrid.jl` is defined by a `Grid` and a `Stratigraphy`, constructed top-down from individual `Layer`s and `Process`es. 

```@docs
variables(::Layer)
variables(::Layer, ::Process)
```
