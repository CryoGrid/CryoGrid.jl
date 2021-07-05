# System type is defined in core/types.jl

# Default empty implementations of diagnostic_step! and prognostic_step! for boundary layers
CryoGrid.diagnosticstep!(::Top, ::System, state) = nothing
CryoGrid.prognosticstep!(::Top, ::System, state) = nothing
CryoGrid.diagnosticstep!(::Bottom, ::System, state) = nothing
CryoGrid.prognosticstep!(::Bottom, ::System, state) = nothing
CryoGrid.variables(layer::Layer, ps::System) = tuplejoin((CryoGrid.variables(layer,p) for p in ps.processes)...)
"""
    CryoGrid.interact!(l1::Layer, ps1::System{P1}, l2::Layer, ps2::System{P2}, s1, s2) where {P1,P2}

Default implementation of `interact!` for multi-process (System) types. Generates a specialized implementation that calls
`interact!` on all pairs of processes between the two layers. Since it is a generated function, the process matching
occurs at compile-time and the emitted code will simply be a sequence of `interact!` calls. Pairs of processes which
lack a definition of `interact!` should be automatically omitted by the compiler.
"""
@generated function CryoGrid.interact!(l1::Layer, ps1::System{P1}, l2::Layer, ps2::System{P2}, s1, s2) where {P1,P2}
    p1types = Tuple(P1.parameters)
    p2types = Tuple(P2.parameters)
    crossprocesses = [(i,j) for (i,p1) in enumerate(p1types) for (j,p2) in enumerate(p2types)]
    expr = Expr(:block)
    for (i,j) in crossprocesses
        @>> quote
        CryoGrid.interact!(l1,ps1[$i],l2,ps2[$j],s1,s2)
        end push!(expr.args)
    end
    return expr
end
"""
    CryoGrid.diagnosticstep!(l::Layer, ps::System{P}, state) where {P}

Default implementation of `diagnosticstep!` for multi-process types. Calls each process in sequence.
"""
@generated function CryoGrid.diagnosticstep!(l::Layer, ps::System{P}, state) where {P}
    expr = Expr(:block)
    for i in 1:length(P.parameters)
        @>> quote
        CryoGrid.diagnosticstep!(l,ps[$i],state)
        end push!(expr.args)
    end
    return expr
end
"""
    CryoGrid.prognosticstep!(l::Layer, ps::System{P}, state) where {P}

Default implementation of `prognosticstep!` for multi-process types. Calls each process in sequence.
"""
@generated function CryoGrid.prognosticstep!(l::Layer, ps::System{P}, state) where {P}
    expr = Expr(:block)
    for i in 1:length(P.parameters)
        @>> quote
        CryoGrid.prognosticstep!(l,ps[$i],state)
        end push!(expr.args)
    end
    return expr
end
"""
    CryoGrid.initialcondition!(l::Layer, ps::System{P}, state) where {P}

Default implementation of `initialcondition!` for multi-process types. Calls each process in sequence.
"""
@generated function CryoGrid.initialcondition!(l::Layer, ps::System{P}, state) where {P}
    expr = Expr(:block)
    for i in 1:length(P.parameters)
        @>> quote
        CryoGrid.initialcondition!(l,ps[$i],state)
        end push!(expr.args)
    end
    return expr
end
"""
    CryoGrid.initialcondition!(l::Layer, ps::System{P}, state) where {P}

Default implementation of `initialcondition!` for multi-process types. Calls each process in sequence.
"""
@generated function CryoGrid.initialcondition!(l1::Layer, ps1::System{P1}, l2::Layer, ps2::System{P2}, s1, s2) where {P1,P2}
    p1types = Tuple(P1.parameters)
    p2types = Tuple(P2.parameters)
    crossprocesses = [(i,j) for (i,p1) in enumerate(p1types) for (j,p2) in enumerate(p2types)]
    expr = Expr(:block)
    for (i,j) in crossprocesses
        @>> quote
        CryoGrid.initialcondition!(l1,ps1[$i],l2,ps2[$j],s1,s2)
        end push!(expr.args)
    end
    return expr
end
