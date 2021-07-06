## Installation

`CryoGrid.jl` can be installed via the Julia package manager:

```
add https://gitlab.awi.de/sparcs/cryogrid/cryogridjulia
```

or equivalently in code/REPL:

```julia
import Pkg
Pkg.add(["https://gitlab.awi.de/sparcs/cryogrid/cryogridjulia"])
```

Be aware that `CryoGrid.jl` is a relatively large package with quite a few dependencies, so installation into a blank Julia environment could take several minutes.

It is recommended that you work with `CryoGrid.jl` as a Julia package rather than cloning the repository and hacking on it directly. This will allow for more rapid development and minimize latency from precompile time. It is also recommended to create a dedicated Julia environment in your workspace to better manage package dependencies. This can be accomplished by running:

```
activate .
```

in your working directory, or by starting Julia with the `--project=.` option. Then, you can proceed to install `CryoGrid.jl` into the environment via the commands above.

You can load `CryoGrid.jl` in your Julia REPL or editor by running:

```julia
using CryoGrid
```

or similarly:

```julia
import CryoGrid
```

The latter option will bring only the `CryoGrid` module name into scope rather than all of its exported components.
