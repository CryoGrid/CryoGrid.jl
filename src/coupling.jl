variables(ps::CoupledProcesses) = tuplejoin((variables(p) for p in ps.process)...)
variables(layer::Layer, ps::CoupledProcesses) = tuplejoin((variables(layer,p) for p in ps.processes)...)
events(layer::Layer, ps::CoupledProcesses) = tuplejoin((events(layer,p) for p in ps.processes)...)
# invoke `f` on all processes in `ps` sequentially
function _invoke_sequential(f!::F, l::Layer, ps::CoupledProcesses, args...) where {F}
    Utils.fastiterate(ps.processes) do proc
        f!(l, proc, args...)
    end
end
CryoGrid.diagnosticstep!(top::Top, ps::CoupledProcesses, state) = _invoke_sequential(diagnosticstep!, top, ps, state)
CryoGrid.diagnosticstep!(bot::Bottom, ps::CoupledProcesses, state) = _invoke_sequential(diagnosticstep!, bot, ps, state)
CryoGrid.diagnosticstep!(sub::SubSurface, ps::CoupledProcesses, state) = _invoke_sequential(diagnosticstep!, sub, ps, state)
CryoGrid.prognosticstep!(top::Top, ps::CoupledProcesses, state) = _invoke_sequential(prognosticstep!, top, ps, state)
CryoGrid.prognosticstep!(bot::Bottom, ps::CoupledProcesses, state) = _invoke_sequential(prognosticstep!, bot, ps, state)
CryoGrid.prognosticstep!(sub::SubSurface, ps::CoupledProcesses, state) = _invoke_sequential(prognosticstep!, sub, ps, state)
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
CryoGrid.initialcondition!(layer::Layer, ps::CoupledProcesses, state) = _invoke_sequential(initialcondition!, layer, ps, state)
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
"""
    timestep(l::Layer, ps::CoupledProcesses{P}, state) where {P}

Default implementation of `timestep` for coupled process types. Calls each process in sequence.
"""
@generated function timestep(l::Layer, ps::CoupledProcesses{P}, state) where {P}
    expr = Expr(:block)
    push!(expr.args, :(dtmax = Inf))
    for i in 1:length(P.parameters)
        quote
            dtmax = min(dtmax, timestep(l, ps[$i], state))
        end |> Base.Fix1(push!, expr.args)
    end
    push!(expr.args, :(return dtmax))
    return expr
end
