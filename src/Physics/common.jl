Constants() = (
    ρw = 1000.0u"kg/m^3", # density of water at standard conditions
    Lsl = 3.34e5u"J/kg", # specific latent heat of fusion of water [J/kg]
    Lsg = 2.257e6u"J/kg", # specific latent heat of vaporization of water [J/kg]
    g = 9.80665u"m/s^2", # gravitational constant
)
# Volume material composition
"""
    volumetricfractions(::SubSurface, ::SubSurfaceProcess, state)
    volumetricfractions(::SubSurface, ::SubSurfaceProcess, state, i)

Get the volumetric fractions of each constituent in the volume (at grid cell `i`, if specificed).
All implementations of `volumetricfractions` are expected to obey a semi-consistent order
in the returned `Tuple` of fractions; the first three consituents should always be `θw,θi,θa`,
i.e. water, ice, and air, followed by any number of additional constituents which may be defined
by the specific layer. There is no feasible way to verify that client code actually obeys this
ordering, so be sure to double check your implementation, otherwise this can cause very subtle bugs!
"""
@inline volumetricfractions(::SubSurface, ::SubSurfaceProcess, state) = ()
@inline volumetricfractions(sub::SubSurface, proc::SubSurfaceProcess, state, i) = volumetricfractions(sub, proc, state)
"""
    partial(f, ::Val{:θw}, sub::SubSurface, proc::SubSurfaceProcess, state, i)

Returns a partially applied function `f` which takes liquid water `θw` as an argument and holds all other
volumetric fractions constant. `f` must be a function of the form `f(::Layer, ::Process, θfracs...)` where
`θfracs` are an arbitrary number of constituent volumetric fractions. Note that this method assumes that
`volumetricfractions` and `f` both obey the implicit ordering convention: `(θw, θi, θa, θfracs...)` where
`θfracs` are zero or more additional constituent fractions. The returned method is a closure which has the
following properties available: `θw, θi, θa, θfracs, θwi` where `θwi` refers to the sum of `θw` and `θi`.
"""
function partial(f::F, ::Val{:θw}, sub::SubSurface, proc::SubSurfaceProcess, state, i) where F
    (_, _, θa, θfracs...) = volumetricfractions(sub, proc, state, i)
    θwi = state.θwi[i]
    return function apply(θw)
        return f(sub, proc, θw, θwi - θw, θa, θfracs...)
    end
end
