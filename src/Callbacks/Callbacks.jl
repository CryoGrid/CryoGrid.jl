module Callbacks

import CryoGrid

using CryoGrid.Drivers: HeatOnlyLandModel
using CryoGrid.Numerics
using CryoGrid.Land: Tile, getvar
using CryoGrid.Utils

using DiffEqCallbacks
using IfElse

"""
CryoGridCallbackFunction{TState,TSetup}(setup, state)

Helper type for defining callbacks on CryoGrid models. Given a Tile and some additional user-defined state type
`TState`, the user can provide dispatches for `CryoGridCallbackFunction{TState}` that satisfy the relevant `DifferentialEquations.jl`
callback function signature. For example:

```julia
struct MyState
    # some state variables
    ...
end
function (fn::CryoGridCallbackFunction{MyState})(u,p,t)
    ...
end
function MyCallback(setup::Tile)
    state = MyState(...)
    fn = CryoGridCallbackFunction(setup, state)
    # create and return SciML callback here
end
```
"""
struct CryoGridCallbackFunction{TState,TSetup}
    setup::TSetup
    state::TState
    CryoGridCallbackFunction(setup::Tile, state::TState) where TState = new{TState, typeof(setup)}(setup, state)
end
(fn::CryoGridCallbackFunction)(u,p,t) = error("no method dispatch provided for callback function $(typeof(fn))")

export CryoGridCallbackFunction

include("courant_step.jl")

end
