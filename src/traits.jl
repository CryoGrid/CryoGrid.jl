# Boundary condition trait
"""
    BoundaryCondition

Trait that specifies the kind of boundary condition. This can be used to write generic
implementations of `interact!` that are (relatively) agnostic to specific implementations of
`BoundaryProcess`. A good example of this can be found in the `boundaryflux` method interface.
"""
abstract type BoundaryCondition end
"""
`BoundaryCondition` instance for Dirichlet boundary conditions.
"""
struct Dirichlet <: BoundaryCondition end
"""
`BoundaryCondition` instance for Neumann boundary conditions.
"""
struct Neumann <: BoundaryCondition end
"""
    BoundaryCondition(::Type{T})

Can be overriden by `BoundaryProcess` types to indicate the type of boundary condition, e.g:

```
BoundaryCondition(::Type{BP}) = Dirichlet()
```

where `BP` is a `BoundaryProcess` that provides the boundary conditions.
"""
BoundaryCondition(::Type{BP}) where {BP<:BoundaryProcess} = error("No boundary condition type specified for boundary process $BP")
BoundaryCondition(bc::BoundaryProcess) = BoundaryCondition(typeof(bc))

"""
    hasfixedvolume(::Type{<:Layer})

Returns true if the given type has fixed/constant volume or false otherwise.
"""
hasfixedvolume(::Type{<:Layer}) = true
