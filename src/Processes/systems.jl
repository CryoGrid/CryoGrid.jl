# Default empty implementations of diagnostic_step! and prognostic_step! for boundary layers
diagnosticstep!(::Top, ::System, state) = nothing
prognosticstep!(::Top, ::System, state) = nothing
diagnosticstep!(::Bottom, ::System, state) = nothing
prognosticstep!(::Bottom, ::System, state) = nothing
variables(layer::Layer, ps::System) = tuplejoin((variables(layer,p) for p in ps.processes)...)
"""
    interact!(l1::Layer, ps1::System{P1}, l2::Layer, ps2::System{P2}, s1, s2) where {P1,P2}

Default implementation of `interact!` for multi-process (System) types. Generates a specialized implementation that calls
`interact!` on all pairs of processes between the two layers. Since it is a generated function, the process matching
occurs at compile-time and the emitted code will simply be a sequence of `interact!` calls. Pairs of processes which
lack a definition of `interact!` should be automatically omitted by the compiler.
"""
@generated function interact!(l1::Layer, ps1::System{P1}, l2::Layer, ps2::System{P2}, s1, s2) where {P1,P2}
    p1types = Tuple(P1.parameters)
    p2types = Tuple(P2.parameters)
    crossprocesses = [(i,j) for (i,p1) in enumerate(p1types) for (j,p2) in enumerate(p2types)]
    expr = Expr(:block)
    for (i,j) in crossprocesses
        @>> quote
        interact!(l1,ps1[$i],l2,ps2[$j],s1,s2)
        end push!(expr.args)
    end
    return expr
end
"""
    diagnosticstep!(l::Layer, ps::System{P}, state) where {P}

Default implementation of `diagnosticstep!` for multi-process types. Calls each process in sequence.
"""
@generated function diagnosticstep!(l::Layer, ps::System{P}, state) where {P}
    expr = Expr(:block)
    for i in 1:length(P.parameters)
        @>> quote
        diagnosticstep!(l,ps[$i],state)
        end push!(expr.args)
    end
    return expr
end
"""
    prognosticstep!(l::Layer, ps::System{P}, state) where {P}

Default implementation of `prognosticstep!` for multi-process types. Calls each process in sequence.
"""
@generated function prognosticstep!(l::Layer, ps::System{P}, state) where {P}
    expr = Expr(:block)
    for i in 1:length(P.parameters)
        @>> quote
        prognosticstep!(l,ps[$i],state)
        end push!(expr.args)
    end
    return expr
end
"""
    initialcondition!(l::Layer, ps::System{P}, state) where {P}

Default implementation of `initialcondition!` for multi-process types. Calls each process in sequence.
"""
@generated function initialcondition!(l::Layer, ps::System{P}, state) where {P}
    expr = Expr(:block)
    for i in 1:length(P.parameters)
        @>> quote
        initialcondition!(l,ps[$i],state)
        end push!(expr.args)
    end
    return expr
end
"""
    initialcondition!(l::Layer, ps::System{P}, state) where {P}

Default implementation of `initialcondition!` for multi-process types. Calls each process in sequence.
"""
@generated function initialcondition!(l1::Layer, ps1::System{P1}, l2::Layer, ps2::System{P2}, s1, s2) where {P1,P2}
    p1types = Tuple(P1.parameters)
    p2types = Tuple(P2.parameters)
    crossprocesses = [(i,j) for (i,p1) in enumerate(p1types) for (j,p2) in enumerate(p2types)]
    expr = Expr(:block)
    for (i,j) in crossprocesses
        @>> quote
        initialcondition!(l1,ps1[$i],l2,ps2[$j],s1,s2)
        end push!(expr.args)
    end
    return expr
end
"""
    observe(::Val{name}, l::Layer, ps::System{P}, state) where {P}

Default implementation of `observe` for multi-process types. Calls each process in sequence.
"""
@generated function observe(val::Val{name}, l::Layer, ps::System{P}, state) where {name,P}
    expr = Expr(:block)
    for i in 1:length(P.parameters)
        @>> quote
        observe(val,l,ps[$i],state)
        end push!(expr.args)
    end
    return expr
end