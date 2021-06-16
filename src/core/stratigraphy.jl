struct StratNode{TLayer,TProcess,name}
    layer::TLayer
    process::TProcess
    StratNode(name::Symbol, layer::TLayer, process::TProcess) where {TLayer<:Layer,TProcess<:Processes} =
        new{TLayer,TProcess,name}(layer,process)
end
"""
Get the name of the given stratigraphy node.
"""
nodename(::StratNode{L,P,name}) where {L,P,name} = name
nodename(::Type{<:StratNode{L,P,name}}) where {L,P,name} = name

Base.show(io::IO, node::StratNode{L,P,name}) where {L,P,name} = print(io, "$name($L,$P)")

# Constructors for stratigraphy nodes
Top(boundaries::BoundaryProcess...) = StratNode(:top, Top(), Processes(boundaries...))
Bottom(boundaries::BoundaryProcess...) = StratNode(:bottom, Bottom(), Processes(boundaries...))
Ground(name::Symbol, layer::SubSurface, processes::SubSurfaceProcess...) = StratNode(name, layer, Processes(processes...))

"""
    Stratigraphy{TNodes,TBounds}

Defines a 1-dimensional stratigraphy by connecting a top and bottom layer to 1 or more subsurface layers.
"""
struct Stratigraphy{TNodes,TBounds}
    boundaries::TBounds
    nodes::TNodes
    Stratigraphy(top::Pair{Q,<:StratNode{Top}}, sub::Pair{Q,<:StratNode{<:SubSurface}},
        bot::Pair{Q,<:StratNode{Bottom}}) where Q = Stratigraphy(top,(sub,),bot)
    Stratigraphy(top::Pair{Q,<:StratNode{Top}}, sub::Tuple{Vararg{Pair{Q,<:StratNode{<:SubSurface}}}},
        bot::Pair{Q,<:StratNode{Bottom}}) where Q = begin
        @assert length(sub) > 0 "At least one subsurface layer must be specified"
        names = @>> sub map(last) map(nodename)
        @assert length(unique(names)) == length(names) "All layer names in Stratigraphy must be unique"
        boundaries = Tuple(map(first, (top, sub..., bot)))
        @assert issorted(boundaries) "Node boundary locations must be in strictly increasing order."
        nodes = Tuple(map(last, (top, sub..., bot)))
        new{typeof(nodes),typeof(boundaries)}(boundaries,nodes)
    end
end
nodetypes(::Type{<:Stratigraphy{TNodes}}) where {TNodes} = Tuple(TNodes.parameters)
# Array and iteration overrides
Base.size(strat::Stratigraphy) = size(strat.nodes)
Base.length(strat::Stratigraphy) = length(strat.nodes)
Base.getindex(strat::Stratigraphy, i::Int) = strat.nodes[i]
Base.iterate(strat::Stratigraphy) = (strat.nodes[1],strat.nodes[2:end])
Base.iterate(strat::Stratigraphy, itrstate::Tuple) = (itrstate[1],itrstate[2:end])
Base.iterate(strat::Stratigraphy, itrstate::Tuple{}) = nothing
Base.show(io::IO, strat::Stratigraphy) = print(io, "Stratigraphy($(prod(("$b => $n, " for (n,b) in zip(strat.nodes,strat.boundaries))))")

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
export Top, Bottom, Ground
export nodename, nodetypes
