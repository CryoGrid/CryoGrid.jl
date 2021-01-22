# Default empty implementations of diagnostic_step! and prognostic_step! for boundary layers
diagnostic_step!(::Top, ::Process, ::State) = nothing
prognostic_step!(::Top, ::Process, ::State) = nothing
diagnostic_step!(::Bottom, ::Process, ::State) = nothing
prognostic_step!(::Bottom, ::Process, ::State) = nothing
@generated function interact!(l1::Layer, ::Processes{P1}, l2::Layer, ::Processes{P2}, s1, s2) where {P1,P2}
    p1types = Tuple(P1.parameters)
    p2types = Tuple(P2.parameters)
    matched_processes = [(p1,p2) for p1 in p1types for p2 in p2types if interactrule(p1,p2)]
    expr = Expr(:block)
    for (p1,p2) in matched_processes
        @>> quote
        interact!(l1,p1,l2,p2,s1,s2)
        end push!(expr.args)
    end
    return expr
end

include("heat_conduction/heat_conduction.jl")
