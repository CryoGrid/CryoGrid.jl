"""
    StratComponent{TLayer,TProcess,name}

Represents a single component (layer + processes) in the stratigraphy.
"""
struct StratComponent{TLayer,TProcess,name}
    layer::TLayer
    process::TProcess
    StratComponent(name::Symbol, layer::TLayer, process::TProcess) where {TLayer<:Layer,TProcess<:CompoundProcess} =
        new{TLayer,TProcess,name}(layer,process)
end
ConstructionBase.constructorof(::Type{StratComponent{TLayer,TProcess,name}}) where {TLayer,TProcess,name} = (layer,process) -> StratComponent(name, layer, process)
"""
Get the name of the given stratigraphy node.
"""
componentname(::StratComponent{L,P,name}) where {L,P,name} = name
componentname(::Type{<:StratComponent{L,P,name}}) where {L,P,name} = name

Base.show(io::IO, node::StratComponent{L,P,name}) where {L,P,name} = print(io, "$name($L,$P)")

# Constructors for stratigraphy nodes
top(bcs::BoundaryProcess...) = StratComponent(:top, Top(), CompoundProcess(bcs...))
bottom(bcs::BoundaryProcess...) = StratComponent(:bottom, Bottom(), CompoundProcess(bcs...))
subsurface(name::Symbol, layer::SubSurface, processes::SubSurfaceProcess...) = StratComponent(name, layer, CompoundProcess(processes...))

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
boundaryintervals(strat::Stratigraphy, z_bottom) = boundaryintervals(boundaries(strat), z_bottom)
boundaryintervals(bounds::NTuple, z_bottom) = tuplejoin(map(tuple, bounds[1:end-1], bounds[2:end]), ((bounds[end], z_bottom),))
componenttypes(::Type{<:Stratigraphy{N,TComponents}}) where {N,TComponents} = Tuple(TComponents.parameters)
Base.getproperty(strat::Stratigraphy, sym::Symbol) = strat[Val{sym}()]
# Array and iteration overrides
Base.size(strat::Stratigraphy) = size(components(strat))
Base.length(strat::Stratigraphy) = length(components(strat))
Base.getindex(strat::Stratigraphy, i::Int) = components(strat)[i]
Base.getindex(strat::Stratigraphy, sym::Symbol) = strat[Val{sym}]
Base.getindex(strat::Stratigraphy, ::Val{sym}) where {sym} = strat[Val{sym}]
@generated Base.getindex(strat::Stratigraphy{N,TC}, ::Type{Val{sym}}) where {N,TC,sym} = :(components(strat)[$(findfirst(T -> componentname(T) == sym, TC.parameters))])
Base.iterate(strat::Stratigraphy) = (components(strat)[1],components(strat)[2:end])
Base.iterate(strat::Stratigraphy, itrstate::Tuple) = (itrstate[1],itrstate[2:end])
Base.iterate(strat::Stratigraphy, itrstate::Tuple{}) = nothing
Base.show(io::IO, strat::Stratigraphy) = print(io, "Stratigraphy($(prod(("$b => $n, " for (n,b) in zip(components(strat),boundaries(strat)))))")

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
