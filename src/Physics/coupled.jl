# Default empty implementations of diagnostic_step! and prognostic_step! for boundary layers
CryoGrid.diagnosticstep!(::Top, ::BoundaryProcess, state) = nothing
CryoGrid.prognosticstep!(::Top, ::BoundaryProcess, state) = nothing
CryoGrid.diagnosticstep!(::Bottom, ::BoundaryProcess, state) = nothing
CryoGrid.prognosticstep!(::Bottom, ::BoundaryProcess, state) = nothing
CryoGrid.variables(ps::CoupledProcesses) = tuplejoin((CryoGrid.variables(p) for p in ps.processes)...)
CryoGrid.variables(layer::Layer, ps::CoupledProcesses) = tuplejoin((CryoGrid.variables(layer,p) for p in ps.processes)...)
CryoGrid.events(layer::Layer, ps::CoupledProcesses) = tuplejoin((CryoGrid.events(layer,p) for p in ps.processes)...)
"""
    interact!(l1::Layer, ps1::CoupledProcesses{P1}, l2::Layer, ps2::CoupledProcesses{P2}, s1, s2) where {P1,P2}

Default implementation of `interact!` for multi-process (CoupledProcesses) types. Generates a specialized implementation that calls
`interact!` on all pairs of processes between the two layers. Since it is a generated function, the process matching
occurs at compile-time and the emitted code will simply be a sequence of `interact!` calls. Pairs of processes which
lack a definition of `interact!` should be automatically omitted by the compiler.
"""
@generated function CryoGrid.interact!(l1::Layer, ps1::CoupledProcesses{P1}, l2::Layer, ps2::CoupledProcesses{P2}, s1, s2) where {P1,P2}
    expr = Expr(:block)
    for i in 1:length(P1.parameters)
        @>> quote
        CryoGrid.interact!(l1,ps1[$i],l2,ps2,s1,s2)
        end push!(expr.args)
    end
    for i in 1:length(P2.parameters)
        @>> quote
        CryoGrid.interact!(l1,ps1,l2,ps2[$i],s1,s2)
        end push!(expr.args)
    end
    return expr
end
@generated function CryoGrid.interact!(l1::Layer, p1::Process, l2::Layer, ps2::CoupledProcesses{P2}, s1, s2) where {P2}
    expr = Expr(:block)
    for i in 1:length(P2.parameters)
        @>> quote
        CryoGrid.interact!(l1,p1,l2,ps2[$i],s1,s2)
        end push!(expr.args)
    end
    return expr
end
@generated function CryoGrid.interact!(l1::Layer, ps1::CoupledProcesses{P1}, l2::Layer, p2::Process, s1, s2) where {P1}
    expr = Expr(:block)
    for i in 1:length(P1.parameters)
        @>> quote
        CryoGrid.interact!(l1,ps1[$i],l2,p2,s1,s2)
        end push!(expr.args)
    end
    return expr
end
"""
    diagnosticstep!(l::Layer, ps::CoupledProcesses{P}, state) where {P}

Default implementation of `diagnosticstep!` for multi-process types. Calls each process in sequence.
"""
@generated function CryoGrid.diagnosticstep!(l::Layer, ps::CoupledProcesses{P}, state) where {P}
    expr = Expr(:block)
    for i in 1:length(P.parameters)
        @>> quote
        CryoGrid.diagnosticstep!(l,ps[$i],state)
        end push!(expr.args)
    end
    return expr
end
"""
    prognosticstep!(l::Layer, ps::CoupledProcesses{P}, state) where {P}

Default implementation of `prognosticstep!` for multi-process types. Calls each process in sequence.
"""
@generated function CryoGrid.prognosticstep!(l::Layer, ps::CoupledProcesses{P}, state) where {P}
    expr = Expr(:block)
    for i in 1:length(P.parameters)
        @>> quote
        CryoGrid.prognosticstep!(l,ps[$i],state)
        end push!(expr.args)
    end
    return expr
end
"""
    initialcondition!([l::Layer,] ps::CoupledProcesses{P}, state) where {P}

Default implementation of `initialcondition!` for multi-process types. Calls each process in sequence.
"""
@generated function CryoGrid.initialcondition!(ps::CoupledProcesses{P}, state) where {P}
    expr = Expr(:block)
    for i in 1:length(P.parameters)
        @>> quote
        CryoGrid.initialcondition!(ps[$i],state)
        end push!(expr.args)
    end
    return expr
end
@generated function CryoGrid.initialcondition!(l::Layer, ps::CoupledProcesses{P}, state) where {P}
    expr = Expr(:block)
    for i in 1:length(P.parameters)
        @>> quote
        CryoGrid.initialcondition!(l,ps[$i],state)
        end push!(expr.args)
    end
    return expr
end
"""
    initialcondition!(l1::Layer, ps1::CoupledProcesses{P1}, l2::Layer, ps2::CoupledProcesses{P2}, s1, s2) where {P1,P2}

Default implementation of `initialcondition!` for multi-process types. Calls each process in sequence.
"""
@generated function CryoGrid.initialcondition!(l1::Layer, ps1::CoupledProcesses{P1}, l2::Layer, ps2::CoupledProcesses{P2}, s1, s2) where {P1,P2}
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
@generated function CryoGrid.criterion(ev::ContinuousEvent, l::Layer, ps::CoupledProcesses{P}, state) where {name,P}
    expr = Expr(:block)
    push!(expr.args, :(value = 1.0))
    for i in 1:length(P.parameters)
        @>> quote
        value *= CryoGrid.criterion(ev,l,ps[$i],state)
        end push!(expr.args)
    end
    push!(expr.args, :(return value))
    return expr
end
@generated function CryoGrid.criterion(ev::DiscreteEvent, l::Layer, ps::CoupledProcesses{P}, state) where {name,P}
    expr = Expr(:block)
    push!(expr.args, :(value = true))
    for i in 1:length(P.parameters)
        @>> quote
        value *= CryoGrid.criterion(ev,l,ps[$i],state)
        end push!(expr.args)
    end
    push!(expr.args, :(return value))
    return expr
end
@generated function CryoGrid.trigger!(ev::Event, l::Layer, ps::CoupledProcesses{P}, state) where {name,P}
    expr = Expr(:block)
    for i in 1:length(P.parameters)
        @>> quote
        CryoGrid.trigger!(ev,l,ps[$i],state)
        end push!(expr.args)
    end
    return expr
end
@generated function CryoGrid.trigger!(ev::ContinuousEvent, tr::ContinuousTrigger, l::Layer, ps::CoupledProcesses{P}, state) where {name,P}
    expr = Expr(:block)
    for i in 1:length(P.parameters)
        @>> quote
        CryoGrid.trigger!(ev,tr,l,ps[$i],state)
        end push!(expr.args)
    end
    return expr
end
"""
    observe(::Val{name}, l::Layer, ps::CoupledProcesses{P}, state) where {P}

Default implementation of `observe` for multi-process types. Calls each process in sequence.
"""
@generated function CryoGrid.observe(val::Val{name}, l::Layer, ps::CoupledProcesses{P}, state) where {name,P}
    expr = Expr(:block)
    for i in 1:length(P.parameters)
        @>> quote
        CryoGrid.observe(val,l,ps[$i],state)
        end push!(expr.args)
    end
    return expr
end