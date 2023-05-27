import Flatten: flattenable

export volumetricfractions

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

export ConstantBC, PeriodicBC
export ConstantValue, PeriodicValue, ConstantFlux, PeriodicFlux
include("simple_bc.jl")
include("composite_bc.jl")
include("steplimiters.jl")
# Sub modules
include("Hydrology/Hydrology.jl")
@reexport using .Hydrology
include("Heat/Heat.jl")
@reexport using .Heat
include("Snow/Snow.jl")
@reexport using .Snow
include("Soils/Soils.jl")
@reexport using .Soils
include("Salt/Salt.jl")
@reexport using .Salt
include("Surface/Surface.jl")
@reexport using .Surface
include("Sources/Sources.jl")
@reexport using .Sources
