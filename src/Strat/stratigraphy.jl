"""
    StratComponent{TLayer,TProcess,name}

Represents a single component (layer + processes) in the stratigraphy.
"""
struct StratComponent{TLayer,TProcesses,name}
    layer::TLayer
    processes::TProcesses
    StratComponent(name::Symbol, layer::TLayer, processes::TProcesses) where {TLayer<:Layer,TProcesses<:CoupledProcesses} =
        new{TLayer,TProcesses,name}(layer,processes)
end
ConstructionBase.constructorof(::Type{StratComponent{TLayer,TProcesses,name}}) where {TLayer,TProcesses,name} = (layer,processes) -> StratComponent(name, layer, processes)
"""
Get the name of the given stratigraphy node.
"""
componentname(::StratComponent{L,P,name}) where {L,P,name} = name
componentname(::Type{<:StratComponent{L,P,name}}) where {L,P,name} = name

Base.show(io::IO, node::StratComponent{L,P,name}) where {L,P,name} = print(io, "$name($L,$P)")

# Constructors for stratigraphy nodes
top(bcs::BoundaryProcess...) = StratComponent(:top, Top(), CoupledProcesses(bcs...))
bottom(bcs::BoundaryProcess...) = StratComponent(:bottom, Bottom(), CoupledProcesses(bcs...))
subsurface(name::Symbol, sub::SubSurface, processes::SubSurfaceProcess...) = StratComponent(name, sub, CoupledProcesses(processes...))

"""
    Stratigraphy{N,TComponents,Q}

Defines a 1-dimensional stratigraphy by connecting a top and bottom layer to 1 or more subsurface layers.
"""
struct Stratigraphy{N,TComponents,Q}
    boundaries::NTuple{N,Q}
    components::TComponents
    Stratigraphy(boundaries::NTuple{N,Q}, components::NTuple{N,StratComponent}) where {N,Q} = new{N,typeof(components),Q}(boundaries, components)
    Stratigraphy(top::Pair{Q,<:StratComponent{Top}}, sub::Pair{Q,<:StratComponent{<:SubSurface}},
        bot::Pair{Q,<:StratComponent{Bottom}}) where {Q<:DistQuantity} = Stratigraphy(top,(sub,),bot)
    function Stratigraphy(
        # use @nospecialize to (hopefully) reduce compilation overhead
        @nospecialize(top::Pair{Q,<:StratComponent{Top}}),
        @nospecialize(sub::Tuple{Vararg{Pair{Q,<:StratComponent{<:SubSurface}}}}),
        @nospecialize(bot::Pair{Q,<:StratComponent{Bottom}})
    ) where {Q<:DistQuantity}
        @assert length(sub) > 0 "At least one subsurface layer must be specified"
        names = map(componentname, map(last, sub))
        @assert length(unique(names)) == length(names) "All layer names in Stratigraphy must be unique"
        boundaries = Tuple(map(pair -> Param(ustrip(first(pair)), units=unit(Q), layer=:strat), (top, sub..., bot)))
        @assert issorted(boundaries, by=p -> p.val) "Stratigraphy boundary locations must be in strictly increasing order."
        components = Tuple(map(last, (top, sub..., bot)))
        new{length(components),typeof(components),eltype(boundaries)}(boundaries,components)
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
