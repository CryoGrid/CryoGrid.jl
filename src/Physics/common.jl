Constants() = (
    ρw = 1000.0u"kg/m^3", # density of water at standard conditions
    Lf = 3.34e5u"J/kg", # specific latent heat of fusion of water [J/kg]
    g = 9.80665u"m/s^2", # gravitational constant
)
# Composition
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
    waterice(::SubSurface, state)
    waterice(sub::SubSurface, ::SubSurfaceProcess, state)
    waterice(sub::SubSurface, ::SubSurfaceProcess, state, i)

Retrieves the total water content (water + ice) for the given layer at grid cell `i`, if specified.
Defaults to retrieving the state variable `θwi` (assuming it exists).
"""
@inline waterice(::SubSurface, state) = state.θwi
@inline waterice(sub::SubSurface, proc::SubSurfaceProcess, state) = waterice(sub, state)
@inline waterice(sub::SubSurface, proc::SubSurfaceProcess, state, i) = Utils.getscalar(waterice(sub, proc, state), i)
"""
    liquidwater(sub::SubSurface, proc::SubSurfaceProcess, state)
    liquidwater(sub::SubSurface, proc::SubSurfaceProcess, state, i)

Retrieves the liquid water content for the given layer.
"""
@inline liquidwater(::SubSurface, ::SubSurfaceProcess, state) = state.θw
@inline liquidwater(sub::SubSurface, proc::SubSurfaceProcess, state, i) = Utils.getscalar(liquidwater(sub, proc, state), i)
"""
    partial(f, ::typeof(liquidwater), sub::SubSurface, proc::SubSurfaceProcess, state, i)

Returns a partially applied function `f` which takes liquid water `θw` as an argument and holds all other
volumetric fractions constant. `f` must be a function of the form `f(::Layer, ::Process, θfracs...)` where
`θfracs` are an arbitrary number of constituent volumetric fractions. Note that this method assumes that
`volumetricfractions` and `f` both obey the implicit ordering convention: `(θw, θi, θa, θfracs...)` where
`θfracs` are zero or more additional constituent fractions. The returned method is a closure which has the
following properties available: `θw, θi, θa, θfracs, θwi` where `θwi` refers to the sum of `θw` and `θi`.
"""
function partial(f::F, ::typeof(liquidwater), sub::SubSurface, proc::SubSurfaceProcess, state, i) where F
    (θw, θi, θa, θfracs...) = volumetricfractions(sub, proc, state, i)
    θwi = θw + θi
    return function apply(θw)
        return f(sub, proc, θw, θwi - θw, θa, θfracs...)
    end
end
