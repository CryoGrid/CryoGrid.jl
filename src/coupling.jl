variables(ps::CoupledProcesses) = tuplejoin((variables(p) for p in ps.process)...)
variables(layer::Layer, ps::CoupledProcesses) = tuplejoin((variables(layer,p) for p in ps.processes)...)
events(layer::Layer, ps::CoupledProcesses) = tuplejoin((events(layer,p) for p in ps.processes)...)
"""
    interact!(l1::Layer, ps1::CoupledProcesses{P1}, l2::Layer, ps2::CoupledProcesses{P2}, s1, s2) where {P1,P2}

Default implementation of `interact!` for coupled process (CoupledProcesses) types. Generates a specialized implementation that calls
`interact!` on all pairs of processes between the two layers. Since it is a generated function, the process matching
occurs at compile-time and the emitted code will simply be a sequence of `interact!` calls. Pairs of processes which
lack a definition of `interact!` should be automatically omitted by the compiler.
"""
@generated function interact!(l1::Layer, ps1::CoupledProcesses{P1}, l2::Layer, ps2::CoupledProcesses{P2}, s1, s2) where {P1,P2}
    p1types = Tuple(P1.parameters)
    p2types = Tuple(P2.parameters)
    crossprocesses = [(i,j) for i in 1:length(p1types) for j in 1:length(p2types)]
    expr = Expr(:block)
    for (i,j) in crossprocesses
        quote
            interact!(l1,ps1[$i],l2,ps2[$j],s1,s2)
        end |> Base.Fix1(push!, expr.args)
    end
    return expr
end
@generated function interact!(l1::Layer, p1::Process, l2::Layer, ps2::CoupledProcesses{P2}, s1, s2) where {P2}
    expr = Expr(:block)
    for i in 1:length(P2.parameters)
        quote
            interact!(l1,p1,l2,ps2[$i],s1,s2)
        end |> Base.Fix1(push!, expr.args)
    end
    return expr
end
@generated function interact!(l1::Layer, ps1::CoupledProcesses{P1}, l2::Layer, p2::Process, s1, s2) where {P1}
    expr = Expr(:block)
    for i in 1:length(P1.parameters)
        quote
            interact!(l1,ps1[$i],l2,p2,s1,s2)
        end |> Base.Fix1(push!, expr.args)
    end
    return expr
end
"""
    diagnosticstep!(l::Layer, ps::CoupledProcesses{P}, state) where {P}

Default implementation of `diagnosticstep!` for coupled process types. Calls each process in sequence.
"""
@generated function diagnosticstep!(l::Layer, ps::CoupledProcesses{P}, state) where {P}
    expr = Expr(:block)
    for i in 1:length(P.parameters)
        quote
            diagnosticstep!(l,ps[$i],state)
        end |> Base.Fix1(push!, expr.args)
    end
    return expr
end
"""
    prognosticstep!(l::Layer, ps::CoupledProcesses{P}, state) where {P}

Default implementation of `prognosticstep!` for coupled process types. Calls each process in sequence.
"""
@generated function prognosticstep!(l::Layer, ps::CoupledProcesses{P}, state) where {P}
    expr = Expr(:block)
    for i in 1:length(P.parameters)
        quote
            prognosticstep!(l,ps[$i],state)
        end |> Base.Fix1(push!, expr.args)
    end
    return expr
end
"""
    initialcondition!([l::Layer,] ps::CoupledProcesses{P}, state) where {P}

Default implementation of `initialcondition!` for coupled process types. Calls each process in sequence.
"""
@generated function initialcondition!(ps::CoupledProcesses{P}, state) where {P}
    expr = Expr(:block)
    for i in 1:length(P.parameters)
        quote
            initialcondition!(ps[$i],state)
        end |> Base.Fix1(push!, expr.args)
    end
    return expr
end
@generated function initialcondition!(l::Layer, ps::CoupledProcesses{P}, state) where {P}
    expr = Expr(:block)
    for i in 1:length(P.parameters)
        quote
            initialcondition!(l,ps[$i],state)
        end |> Base.Fix1(push!, expr.args)
    end
    return expr
end
"""
    initialcondition!(l1::Layer, ps1::CoupledProcesses{P1}, l2::Layer, ps2::CoupledProcesses{P2}, s1, s2) where {P1,P2}

Default implementation of `initialcondition!` for coupled process types. Calls each process in sequence.
"""
@generated function initialcondition!(l1::Layer, ps1::CoupledProcesses{P1}, l2::Layer, ps2::CoupledProcesses{P2}, s1, s2) where {P1,P2}
    p1types = Tuple(P1.parameters)
    p2types = Tuple(P2.parameters)
    crossprocesses = [(i,j) for i in 1:length(p1types) for j in 1:length(p2types)]
    expr = Expr(:block)
    for (i,j) in crossprocesses
        quote
            initialcondition!(l1,ps1[$i],l2,ps2[$j],s1,s2)
        end |> Base.Fix1(push!, expr.args)
    end
    return expr
end
@generated function initialcondition!(l1::Layer, p1::Process, l2::Layer, ps2::CoupledProcesses{P2}, s1, s2) where {P2}
    expr = Expr(:block)
    for i in 1:length(P2.parameters)
        quote
            initialcondition!(l1,p1,l2,ps2[$i],s1,s2)
        end |> Base.Fix1(push!, expr.args)
    end
    return expr
end
@generated function initialcondition!(l1::Layer, ps1::CoupledProcesses{P1}, l2::Layer, p2::Process, s1, s2) where {P1}
    expr = Expr(:block)
    for i in 1:length(P1.parameters)
        quote
            initialcondition!(l1,ps1[$i],l2,p2,s1,s2)
        end |> Base.Fix1(push!, expr.args)
    end
    return expr
end
@generated function criterion(ev::ContinuousEvent, l::Layer, ps::CoupledProcesses{P}, state) where {name,P}
    expr = Expr(:block)
    push!(expr.args, :(value = 1.0))
    for i in 1:length(P.parameters)
        quote
            value *= criterion(ev,l,ps[$i],state)
        end |> Base.Fix1(push!, expr.args)
    end
    push!(expr.args, :(return value))
    return expr
end
@generated function criterion(ev::DiscreteEvent, l::Layer, ps::CoupledProcesses{P}, state) where {name,P}
    expr = Expr(:block)
    push!(expr.args, :(value = true))
    for i in 1:length(P.parameters)
        quote
            value *= criterion(ev,l,ps[$i],state)
        end |> Base.Fix1(push!, expr.args)
    end
    push!(expr.args, :(return value))
    return expr
end
@generated function trigger!(ev::Event, l::Layer, ps::CoupledProcesses{P}, state) where {name,P}
    expr = Expr(:block)
    for i in 1:length(P.parameters)
        quote
            trigger!(ev,l,ps[$i],state)
        end |> Base.Fix1(push!, expr.args)
    end
    return expr
end
@generated function trigger!(ev::ContinuousEvent, tr::ContinuousTrigger, l::Layer, ps::CoupledProcesses{P}, state) where {name,P}
    expr = Expr(:block)
    for i in 1:length(P.parameters)
        quote
            trigger!(ev,tr,l,ps[$i],state)
        end |> Base.Fix1(push!, expr.args)
    end
    return expr
end
"""
    timestep(l::Layer, ps::CoupledProcesses{P}, state) where {P}

Default implementation of `timestep` for coupled process types. Calls each process in sequence.
"""
@generated function timestep(l::Layer, ps::CoupledProcesses{P}, state) where {P}
    expr = Expr(:block)
    push!(expr.args, :(dtmax = Inf))
    for i in 1:length(P.parameters)
        quote
            dtmax = min(dtmax, timestep(l,ps[$i],state))
        end |> Base.Fix1(push!, expr.args)
    end
    push!(expr.args, :(return dtmax))
    return expr
end
"""
    observe(::Val{name}, l::Layer, ps::CoupledProcesses{P}, state) where {P}

Default implementation of `observe` for coupled process types. Calls each process in sequence.
"""
@generated function observe(val::Val{name}, l::Layer, ps::CoupledProcesses{P}, state) where {name,P}
    expr = Expr(:block)
    for i in 1:length(P.parameters)
        quote
            observe(val,l,ps[$i],state)
        end |> Base.Fix1(push!, expr.args)
    end
    return expr
end