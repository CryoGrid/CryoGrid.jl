struct StratNode{TLayer,TProcess,name}
    layer::TLayer
    process::TProcess
    StratNode(name::Symbol, layer::TLayer, process::TProcess) where {TLayer<:Layer,TProcess<:Process} =
        new{TLayer,TProcess,name}(layer,process)
end
"""
Get the name of the given stratigraphy node.
"""
name(::StratNode{L,P,name}) where {L,P,I,name} = name
name(::Type{StratNode{L,P,name}}) where {L,P,I,name} = name

# Constructors for stratigraphy nodes
Top(boundaries::BoundaryProcess...) = StratNode(:top, Top(), Processes(boundaries...))
Bottom(boundaries::BoundaryProcess...) = StratNode(:bottom, Bottom(), Processes(boundaries...))
Ground(layer::SubSurface, name::Symbol, processes::SubSurfaceProcess...) = StratNode(name, layer, Processes(processes...))

"""
    Stratigraphy{TNodes,N,Q}

Defines a 1-dimensional stratigraphy by connecting a top and bottom layer to 1 or more subsurface layers.
"""
struct Stratigraphy{TNodes,TBounds}
    boundaries::TBounds
    nodes::TNodes
    Stratigraphy(top::StratNode{<:Top,<:Process},
        sub::Pair{Q,StratNode{<:SubSurface,<:Process}},
        bot::StratNode{<:Bottom,<:Process}) where {Q<:DistQuantity} = Stratigraphy(top,(sub,),bot)
    Stratigraphy(top::StratNode{<:Top,<:Process},
        sub::Tuple{Vararg{Pair{Q,StratNode{<:SubSurface,<:Process}}}},
        bot::Pair{Q,StratNode{<:Bottom,<:Process}}) where {Q<:DistQuantity} = begin
        @assert length(sub) > 0 "At least one subsurface layer must be specified"
        dict = SortedDict(sub..., bot)
        boundaries = tuple(keys(dict)...)
        nodes = tuple(top,values(dict)...)
        new{typeof(nodes),typeof(boundaries)}(boundaries,nodes)
    end
end
nodetypes(::Type{Stratigraphy{TNodes}}) where {TNodes} = Tuple(TNodes.parameters)
# Array and iteration overrides
Base.size(strat::Stratigraphy) = size(strat.nodes)
Base.length(strat::Stratigraphy) = length(strat.nodes)
Base.getindex(strat::Stratigraphy, i::Int) = strat.nodes[i]
Base.iterate(strat::Stratigraphy) = (strat.nodes[1],strat.nodes)
Base.iterate(strat::Stratigraphy, itrstate::Tuple) = (itrstate[1],itrstate[2:end])
Base.iterate(strat::Stratigraphy, itrstate::Tuple{}) = nothing

export Stratigraphy
export Top, Bottom, Ground
export layerat
