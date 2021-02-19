# Default empty implementations of diagnostic_step! and prognostic_step! for boundary layers
diagnostic_step!(::Top, ::Process, ::State) = nothing
prognostic_step!(::Top, ::Process, ::State) = nothing
diagnostic_step!(::Bottom, ::Process, ::State) = nothing
prognostic_step!(::Bottom, ::Process, ::State) = nothing
"""
Default implementation of `interact!` for multi-process types. Generates a specialized implementation that calls
`interact!` on all compatible pairs of processes between the two layers. Since it is a generated function, the
process matching occurs at compile-time and the emitted code will simply be a sequence of `interact!` calls.
"""
@generated function interact!(l1::Layer, ps1::Processes{P1}, l2::Layer, ps2::Processes{P2}, s1, s2) where {P1,P2}
    p1types = Tuple(P1.parameters)
    p2types = Tuple(P2.parameters)
    matched_processes = [(i,j) for (i,p1) in enumerate(p1types) for (j,p2) in enumerate(p2types) if interactrule(p1,p2)]
    expr = Expr(:block)
    for (i,j) in matched_processes
        @>> quote
        interact!(l1,ps1[$i],l2,ps2[$j],s1,s2)
        end push!(expr.args)
    end
    return expr
end

include("boundary/boundaries.jl")
include("heat_conduction/heat_conduction.jl")
