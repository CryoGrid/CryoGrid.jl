module Physics

import CryoGrid

using CryoGrid
using CryoGrid.Numerics
using CryoGrid.Utils

using IfElse
using Reexport
using Unitful

export volumetricfractions, partial

Constants = (
    ρw = 1000.0u"kg/m^3", # density of water at standard conditions
    Lsl = 3.34e5u"J/kg", # specific latent heat of fusion of water [J/kg]
    Lsg = 2.257e6u"J/kg", # specific latent heat of vaporization of water [J/kg]
    g = 9.80665u"m/s^2", # gravitational constant
)

# Volume material composition
"""
    volumetricfractions(::SubSurface, state)
    volumetricfractions(::SubSurface, state, i)

Get the volumetric fractions of each constituent in the volume (at grid cell `i`, if specificed).
All implementations of `volumetricfractions` are expected to obey a semi-consistent order
in the returned `Tuple` of fractions; the first three consituents should always be `θw,θi,θa`,
i.e. water, ice, and air, followed by any number of additional constituents which may be defined
by the specific layer. There is no feasible way to verify that client code actually obeys this
ordering, so be sure to double check your implementation, otherwise this can cause very subtle bugs!
"""
@inline volumetricfractions(::SubSurface, state) = ()
@inline volumetricfractions(sub::SubSurface, state, i) = volumetricfractions(sub, state)
"""
    partial(f, ::Val{:θw}, sub::SubSurface, proc::Process, state, i)

Returns a partially applied function `f` which takes liquid water `θw` as an argument and holds all other
volumetric fractions constant. `f` must be a function of the form `f(::Layer, ::Process, θfracs...)` where
`θfracs` are an arbitrary number of constituent volumetric fractions. Note that this method assumes that
`volumetricfractions` and `f` both obey the implicit ordering convention: `(θw, θi, θa, θfracs...)` where
`θfracs` are zero or more additional constituent fractions. The returned method is a closure which has the
following properties available: `θw, θi, θa, θfracs, θwi` where `θwi` refers to the sum of `θw` and `θi`.
"""
function partial(f::F, ::Val{:θw}, sub::SubSurface, proc::Process, state, i) where F
    function apply(θw)
        (_, _, θa, θfracs...) = volumetricfractions(sub, state, i)
        θwi = state.θwi[i]
        return f(sub, proc, θw, θwi - θw, θa, θfracs...)
    end
end

include("steplimiters.jl")
include("Boundaries/Boundaries.jl")
include("Hydrology/Hydrology.jl")
include("Heat/Heat.jl")
include("Snow/Snow.jl")
include("Soils/Soils.jl")
include("SEB/SEB.jl")
include("Sources/Sources.jl")

@reexport using .Boundaries
@reexport using .Heat
@reexport using .Hydrology
@reexport using .Snow
@reexport using .Soils
@reexport using .SEB
@reexport using .Sources

end