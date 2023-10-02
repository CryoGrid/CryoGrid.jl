# [Concepts](@id concepts)

This page provides a brief(-ish) overview of the basic programming and Julia-language concepts needed in order to develop and extend CryoGrid.jl.
 
## Multiple dispatch

The [CryoGrid community model](https://github.com/CryoGrid/CryoGridCommunity_source) in MATLAB reuses code through [object-oriented programming]() (OOP), namely by separating process implementations into standalone classes. These classes are then subtyped by layer classes (e.g. `GROUND_freeW_ubT`) which then inherit their state variables, parameters, and functions.

Unlike MATLAB, Julia is not object-oriented and has no concept of a class (the closest equivalent is a [struct](https://docs.julialang.org/en/v1/manual/types/#Composite-Types)). Furthermore, while Julia allows for inheritance from [abstract types](https://docs.julialang.org/en/v1/manual/types/#man-abstract-types), it does not allow `struct`s to inherit from other `struct`s, as explained by the [Julia documentation](https://docs.julialang.org/en/v1/manual/types/):

> One particularly distinctive feature of Julia's type system is that concrete types may not subtype each other: all concrete types are final and may only have abstract types as their supertypes. While this might at first seem unduly restrictive, it has many beneficial consequences with surprisingly few drawbacks. It turns out that being able to inherit behavior is much more important than being able to inherit structure, and inheriting both causes significant difficulties in traditional object-oriented languages.

In light of this, CryoGrid.jl takes a different approach to effectively reuse code between model components. Rather than having layers "inherit" processes, subtypes of `Layer` are *composed* of one or more processes as well as parameter types which then determine which methods are invoked at runtime. This is facilitated by one of Julia's key features, [multiple dispatch](https://docs.julialang.org/en/v1/manual/methods/#Defining-Methods). Multiple dispatch means that methods are dynamically invoked based on the (runtime) types of *all* their arguments. This is in constrast to most OOP languages (also MATLAB) where dynamic dispatch occurs based on only one (implicit) argument, i.e. the type of the "object" or class itself. As an example, consider the `CryoGrid` method `updatestate!`:

```julia
using CryoGrid

# Declare two new SubSurface layer types;
# note that <: basically means "is a subtype of"
struct Foo <: SubSurface end
struct Bar <: SubSurface end

CryoGrid.updatestate!(layer::Foo, state) = println("hello Foo")
CryoGrid.updatestate!(layer::Bar, state) = println("hello Bar")

state = nothing # we can ignore the state for the sake of the example
updatestate!(Foo(), state)
updatestate!(Bar(), state)
```
Output:
```
hello Foo
hello Bar
```

In this example, we can see that which `updatestate!` implementation gets invoked is determined by which type is supplied by the caller.

Multiple dispatch allows us to extend this naturally to cases where more than one method argument has a declared type:

```julia
struct MyProcess <: SubSurfaceProcess end

CryoGrid.updatestate!(::SubSurface, ::MyProcess, state) = println("hello MyProcess on any SubSurface")
CryoGrid.updatestate!(::Bar, ::MyProcess, state) = println("hello MyProcess on Bar")

updatestate!(Layer1(), MyProcess() state)
updatestate!(Layer2(), MyProcess(), state)
```
Output:
```
hello MyProcess on any SubSurface
hello MyProcess on Bar
```

Thus, multiple dispatch allows us to write generic code in `updatestate!` that implements `MyProcess` for any `SubSurface` layer (i.e. the parent type of both `Foo` and `Bar`) in addition to adding specialized code for more specific layer types.