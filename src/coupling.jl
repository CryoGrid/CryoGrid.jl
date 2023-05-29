# default variables impl for coupled processes
variables(ps::CoupledProcesses) = tuplejoin((variables(p) for p in ps.process)...)
variables(layer::Layer, ps::CoupledProcesses) = tuplejoin((variables(layer,p) for p in ps.processes)...)

events(layer::Layer, ps::CoupledProcesses) = tuplejoin((events(layer,p) for p in ps.processes)...)

# invoke `f!` on all processes in `ps` sequentially
function _invoke_sequential(f!::F, l::Layer, ps::CoupledProcesses, args...) where {F}
    Utils.fastiterate(ps.processes) do proc
        f!(l, proc, args...)
    end
end

CryoGrid.updatestate!(top::Top, ps::CoupledProcesses, state) = _invoke_sequential(updatestate!, top, ps, state)
CryoGrid.updatestate!(bot::Bottom, ps::CoupledProcesses, state) = _invoke_sequential(updatestate!, bot, ps, state)
CryoGrid.updatestate!(sub::SubSurface, ps::CoupledProcesses, state) = _invoke_sequential(updatestate!, sub, ps, state)

CryoGrid.computefluxes!(top::Top, ps::CoupledProcesses, state) = _invoke_sequential(computefluxes!, top, ps, state)
CryoGrid.computefluxes!(bot::Bottom, ps::CoupledProcesses, state) = _invoke_sequential(computefluxes!, bot, ps, state)
CryoGrid.computefluxes!(sub::SubSurface, ps::CoupledProcesses, state) = _invoke_sequential(computefluxes!, sub, ps, state)

"""
    interact!(l1::Top, ps1::CoupledProcesses{P1}, l2::Layer, ps2::CoupledProcesses{P2}, s1, s2) where {P1,P2}

Default implementation of `interact!` for coupled process (CoupledProcesses) types on the `Top` layer. Iterates over each
boundary process and calls `interact!` for this process and the subsurface layer.
"""
function interact!(l1::Top, ps1::CoupledProcesses, l2::Layer, p2::Process, s1, s2)
    Utils.fastinvoke(ps1) do p_i
        interact!(l1, p_i, l2, p2, s1, s2)
    end
end

"""
    interact!(l1::Layer, p1::Process, l2::Bottom, ps2::CoupledProcesses, s1, s2)

Default implementation of `interact!` for coupled process (CoupledProcesses) types on the `Bottom` layer. Iterates over each
boundary process and calls `interact!` for this process and the subsurface layer.
"""
@generated function interact!(l1::Layer, p1::Process, l2::Bottom, ps2::CoupledProcesses, s1, s2)
    Utils.fastinvoke(ps2) do p_i
        interact!(l1, p1, l2, p_i, s1, s2)
    end
end


"""
    timestep(l::Layer, ps::CoupledProcesses{P}, state) where {P}

Default implementation of `timestep` for coupled process types. Calls each process in sequence.
"""
function timestep(l::Layer, ps::CoupledProcesses{P}, state) where {P}
    dtmax = Utils.fastmap(ps) do p
        timestep(l, p, state)
    end
    return minimum(dtmax)
end
