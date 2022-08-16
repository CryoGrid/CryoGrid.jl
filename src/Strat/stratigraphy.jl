const RESERVED_COMPONENT_NAMES = (:top, :bottom, :strat, :init, :event)

"""
    StratComponent{TLayer,TProc,name}

Represents a single component (layer + process(es)) in the stratigraphy.
"""
struct StratComponent{TLayer,TProc,name}
    layer::TLayer
    process::TProc
    StratComponent(name::Symbol, layer::TLayer, proc::TProc) where {TLayer<:Layer,TProc<:Process} = new{TLayer,TProc,name}(layer, proc)
    function StratComponent(name::Symbol, layer::TLayer, procs::CoupledProcesses; ignore_order=false) where {TLayer<:Layer}
        # check coupling order
        if !issorted(procs.processes) && !ignore_order
            procs = CoupledProcesses(sort(procs.processes))
            @warn "The ordering of the given coupled processes is inconsistent with the defined rules and has been automatically corrected: $(map(p -> typeof(p).name.wrapper, procs.processes)).
            If this was on purpose, you can override the defined ordering and suppress this warning with `ignore_order=true` or define `isless` on your process types to explicitly declare the intended ordering."
        end
        new{TLayer,typeof(procs),name}(layer,procs)
    end
end
ConstructionBase.constructorof(::Type{StratComponent{TLayer,TProc,name}}) where {TLayer,TProc,name} = (layer,process) -> StratComponent(name, layer, process)
"""
Get the name of the given stratigraphy node.
"""
componentname(::StratComponent{L,P,name}) where {L,P,name} = name
componentname(::Type{<:StratComponent{L,P,name}}) where {L,P,name} = name
componentnameval(::StratComponent{L,P,name}) where {L,P,name} = Val{name}

Base.show(io::IO, node::StratComponent{L,P,name}) where {L,P,name} = print(io, "$name($L,$P)")

# Constructors for stratigraphy nodes
top(bc::BoundaryProcess) = StratComponent(:top, Top(), bc)
top(bcs::BoundaryProcess...; ignore_order=false) = StratComponent(:top, Top(), CoupledProcesses(bcs...); ignore_order)
bottom(bc::BoundaryProcess) = StratComponent(:bottom, Bottom(), bc)
bottom(bcs::BoundaryProcess...; ignore_order=false) = StratComponent(:bottom, Bottom(), CoupledProcesses(bcs...); ignore_order)
function subsurface(name::Symbol, sub::SubSurface, proc::SubSurfaceProcess)
    @assert name ∉ RESERVED_COMPONENT_NAMES "layer identifier $name is reserved"
    return StratComponent(name, sub, proc)
end
function subsurface(name::Symbol, sub::SubSurface, procs::SubSurfaceProcess...; ignore_order=false)
    @assert name ∉ RESERVED_COMPONENT_NAMES "layer identifier $name is reserved"
    return StratComponent(name, sub, CoupledProcesses(procs...); ignore_order)
end

"""
Type bound for stratigraphy boundaries. Boundaries may be specified as fixed distance quantities
or via some arbitrary parameterization.
"""
const StratBoundaryType = Union{<:DistQuantity,<:AbstractParam,<:Parameterization}

"""
    Stratigraphy{N,TComponents,TBoundaries}

Defines a 1-dimensional stratigraphy by connecting a top and bottom layer to 1 or more subsurface layers.
"""
struct Stratigraphy{N,TComponents,TBoundaries}
    boundaries::TBoundaries
    components::TComponents
    Stratigraphy(boundaries::NTuple{N,Any}, components::NTuple{N,StratComponent}) where {N} = new{N,typeof(components),typeof(boundaries)}(boundaries, components)
    Stratigraphy(
        top::Pair{<:StratBoundaryType,<:StratComponent{Top}},
        sub::Pair{<:StratBoundaryType,<:StratComponent{<:SubSurface}},
        bot::Pair{<:StratBoundaryType,<:StratComponent{Bottom}}
    ) = Stratigraphy(top,(sub,),bot)
    function Stratigraphy(
        # use @nospecialize to (hopefully) reduce compilation overhead
        @nospecialize(top::Pair{<:StratBoundaryType,<:StratComponent{Top}}),
        @nospecialize(sub::Tuple{Vararg{Pair{<:StratBoundaryType,<:StratComponent{<:SubSurface}}}}),
        @nospecialize(bot::Pair{<:StratBoundaryType,<:StratComponent{Bottom}})
    )
        @assert length(sub) > 0 "At least one subsurface layer must be specified"
        names = map(componentname, map(last, sub))
        @assert length(unique(names)) == length(names) "All layer names in Stratigraphy must be unique"
        boundaries = Tuple(map(first, (top, sub..., bot)))
        @assert issorted(boundaries) "Stratigraphy boundary locations must be in strictly increasing order."
        # get components
        components = Tuple(map(last, (top, sub..., bot)))
        # construct type
        new{length(components),typeof(components),typeof(boundaries)}(boundaries, components)
    end
end
"""
Convenience macro for defining stratigraphies with multiple subsurface layers.
"""
macro Stratigraphy(args...)
    @assert length(args) >= 3 "At least three stratigraphy nodes (top, subsurface, bottom) must be provided!"
    if length(args) == 3
        :(Stratigraphy($(esc(args[1])), $(esc(args[2])), $(esc(args[3]))))
    elseif length(args) > 3
        :(Stratigraphy($(esc(args[1])), tuple($(esc.(args[2:end-1])...)), $(esc(args[end]))))
    end
end
components(strat::Stratigraphy) = getfield(strat, :components)
boundaries(strat::Stratigraphy) = getfield(strat, :boundaries)
boundarypairs(strat::Stratigraphy, z_bottom) = boundarypairs(boundaries(strat), z_bottom)
boundarypairs(bounds::NTuple, z_bottom) = tuplejoin(map(tuple, bounds[1:end-1], bounds[2:end]), ((bounds[end], z_bottom),))
componentnames(strat::Stratigraphy) = map(componentname, components(strat))
componenttypes(::Type{<:Stratigraphy{N,TComponents}}) where {N,TComponents} = Tuple(TComponents.parameters)
Base.keys(strat::Stratigraphy) = componentnames(strat)
Base.values(strat::Stratigraphy) = components(strat)
@inline Base.propertynames(strat::Stratigraphy) = Base.keys(strat)
@inline Base.getproperty(strat::Stratigraphy, sym::Symbol) = strat[Val{sym}()]
@inline Base.getindex(strat::Stratigraphy, sym::Symbol) = strat[Val{sym}()]
@generated Base.getindex(strat::Stratigraphy{N,TC}, ::Val{sym}) where {N,TC,sym} = :(components(strat)[$(findfirst(T -> componentname(T) == sym, TC.parameters))])
# Array and iteration overrides
Base.size(strat::Stratigraphy) = size(components(strat))
Base.length(strat::Stratigraphy) = length(components(strat))
Base.getindex(strat::Stratigraphy, i::Int) = components(strat)[i]
Base.iterate(strat::Stratigraphy) = (components(strat)[1],components(strat)[2:end])
Base.iterate(strat::Stratigraphy, itrstate::Tuple) = (itrstate[1],itrstate[2:end])
Base.iterate(strat::Stratigraphy, itrstate::Tuple{}) = nothing
Base.show(io::IO, strat::Stratigraphy) = print(io, "Stratigraphy($(prod(("$b => $n, " for (n,b) in zip(components(strat),boundaries(strat)))))")
# ConstructionBase
ConstructionBase.getproperties(strat::Stratigraphy) = (;map(Pair, Base.keys(strat), Base.values(strat))...)
function ConstructionBase.setproperties(strat::Stratigraphy, patch::NamedTuple)
    components_patched = map(components(strat)) do comp
        get(patch, componentname(comp), comp)
    end
    return Stratigraphy(boundaries(strat), components_patched)
end
