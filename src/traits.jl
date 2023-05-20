# Boundary condition trait
"""
    BCKind

Trait that specifies the kind of boundary condition. This can be used to write generic
implementations of `interact!` that are (relatively) agnostic to specific implementations of
`BoundaryProcess`. A good example of this can be found in the `boundaryflux` method interface.
"""
abstract type BCKind end
"""
`BCKind` instance for Dirichlet boundary conditions.
"""
struct Dirichlet <: BCKind end
"""
`BCKind` instance for Neumann boundary conditions.
"""
struct Neumann <: BCKind end
"""
    BCKind(::Type{T})

Can be overriden by `BoundaryProcess` types to indicate the type of boundary condition, e.g:

```
BCKind(::Type{BP}) = Dirichlet()
```

where `BP` is a `BoundaryProcess` that provides the boundary conditions.
"""
BCKind(::Type{BP}) where {BP<:BoundaryProcess} = error("No boundary condition type specified for boundary process $BP")
BCKind(bc::BoundaryProcess) = BCKind(typeof(bc))

abstract type Volume end
"""
    FixedVolume

`Volume` trait instance for fixed volume layers. Default for all `Layer` types.
"""
struct FixedVolume <: Volume end
"""
    PrognosticVolume

`Volume` trait instance for layers with time varying volume where the volume should
be treated as a prognostic state variable.
"""
struct PrognosticVolume <: Volume end
"""
    DiagnosticVolume

`Volume` trait instance for layers with time varying volume where the volume should
be treated as a diagnostic state variable.
"""
struct DiagnosticVolume <: Volume end

"""
    Volume(::Type{<:Layer})

Trait for layer types that determines whether its spatial volume is temporally invariant,
`FixedVolume`, or varying with time, `DiagnosticVolume` or `PrognosticVolume`.
"""
Volume(::Type{<:Layer}) = FixedVolume()
Volume(layer::Layer) = Volume(typeof(layer))
