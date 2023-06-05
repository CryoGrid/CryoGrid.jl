# default variables impl for coupled processes
variables(ps::CoupledProcesses) = tuplejoin((variables(p) for p in ps.process)...)
variables(layer::Layer, ps::CoupledProcesses) = tuplejoin((variables(layer,p) for p in ps.processes)...)

events(layer::Layer, ps::CoupledProcesses) = tuplejoin((events(layer,p) for p in ps.processes)...)

initializers(layer::Layer, ps::CoupledProcesses) = tuplejoin((initializers(layer,p) for p in ps.processes)...)

# invoke `f!` on all processes in `ps` sequentially
function _invoke_sequential(f!::F, l::Layer, ps::CoupledProcesses, args...) where {F}
    Utils.fastiterate(ps.processes) do proc
        f!(l, proc, args...)
    end
end

CryoGrid.initialcondition!(layer::Layer, ps::CoupledProcesses, state) = _invoke_sequential(initialcondition!, layer, ps, state)

CryoGrid.updatestate!(layer::Layer, ps::CoupledProcesses, state) = _invoke_sequential(updatestate!, layer, ps, state)

CryoGrid.computefluxes!(layer::Layer, ps::CoupledProcesses, state) = _invoke_sequential(computefluxes!, layer, ps, state)

CryoGrid.resetfluxes!(layer::Layer, ps::CoupledProcesses, state) = _invoke_sequential(resetfluxes!, layer, ps, state)

"""
    timestep(l::Layer, ps::CoupledProcesses{P}, state) where {P}

Default implementation of `timestep` for coupled process types. Calls each process in sequence.
"""
function timestep(l::Layer, ps::CoupledProcesses{P}, state) where {P}
    dtmax = Utils.fastmap(ps.processes) do p
        timestep(l, p, state)
    end
    return minimum(dtmax)
end
