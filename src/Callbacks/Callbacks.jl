module Callbacks

import CryoGrid

using CryoGrid.Numerics
using CryoGrid.Setup: CryoGridSetup, HeatOnlySetup, getvar
using CryoGrid.Utils

using DiffEqCallbacks

"""
CryoGridCallbackFunction{TState,TSetup}(setup, state)

Helper type for defining callbacks on CryoGrid models. Given a CryoGridSetup and some additional user-defined state type
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
function MyCallback(setup::CryoGridSetup)
    state = MyState(...)
    fn = CryoGridCallbackFunction(setup, state)
    # create and return SciML callback here
end
```
"""
struct CryoGridCallbackFunction{TState,TSetup}
    setup::TSetup
    state::TState
    CryoGridCallbackFunction(setup::CryoGridSetup, state::TState) where TState = new{TState, typeof(setup)}(setup, state)
end
(fn::CryoGridCallbackFunction)(u,p,t) = error("no method dispatch provided for callback function $(typeof(fn))")

export CryoGridCallbackFunction

include("courant_step.jl")

end
