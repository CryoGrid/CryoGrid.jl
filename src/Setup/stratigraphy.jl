struct StratComponent{TLayer,TProcess,name}
    layer::TLayer
    process::TProcess
    StratComponent(name::Symbol, layer::TLayer, process::TProcess) where {TLayer<:Layer,TProcess<:CompositeProcess} =
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
top(boundaries::BoundaryProcess...) = StratComponent(:top, Top(), CompositeProcess(boundaries...))
bottom(boundaries::BoundaryProcess...) = StratComponent(:bottom, Bottom(), CompositeProcess(boundaries...))
subsurface(name::Symbol, layer::SubSurface, processes::SubSurfaceProcess...) = StratComponent(name, layer, CompositeProcess(processes...))

"""
    Stratigraphy{N,TComps,TBounds}

Defines a 1-dimensional stratigraphy by connecting a top and bottom layer to 1 or more subsurface layers.
"""
struct Stratigraphy{N,TComps,TBounds}
    boundaries::TBounds
    components::TComps
    Stratigraphy(top::Pair{Q,<:StratComponent{Top}}, sub::Pair{Q,<:StratComponent{<:SubSurface}},
        bot::Pair{Q,<:StratComponent{Bottom}}) where Q = Stratigraphy(top,(sub,),bot)
    Stratigraphy(top::Pair{Q,<:StratComponent{Top}}, sub::Tuple{Vararg{Pair{Q,<:StratComponent{<:SubSurface}}}},
        bot::Pair{Q,<:StratComponent{Bottom}}) where Q = begin
        @assert length(sub) > 0 "At least one subsurface layer must be specified"
        names = @>> sub map(last) map(componentname)
        @assert length(unique(names)) == length(names) "All layer names in Stratigraphy must be unique"
        boundaries = map(first, (top, sub..., bot)) |> Tuple
        @assert issorted(boundaries) "Stratigraphy boundary locations must be in strictly increasing order."
        components = map(last, (top, sub..., bot)) |> Tuple
        new{length(components),typeof(components),typeof(boundaries)}(boundaries,components)
    end
end
copmonenttypes(::Type{<:Stratigraphy{N,TComps}}) where {N,TComps} = Tuple(TComps.parameters)
# Array and iteration overrides
Base.size(strat::Stratigraphy) = size(strat.components)
Base.length(strat::Stratigraphy) = length(strat.components)
Base.getindex(strat::Stratigraphy, i::Int) = strat.components[i]
Base.iterate(strat::Stratigraphy) = (strat.components[1],strat.components[2:end])
Base.iterate(strat::Stratigraphy, itrstate::Tuple) = (itrstate[1],itrstate[2:end])
Base.iterate(strat::Stratigraphy, itrstate::Tuple{}) = nothing
Base.show(io::IO, strat::Stratigraphy) = print(io, "Stratigraphy($(prod(("$b => $n, " for (n,b) in zip(strat.components,strat.boundaries))))")

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

export Stratigraphy, @Stratigraphy
export top, bottom, subsurface
export StratComponent, componentname, copmonenttypes
