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

CryoGrid.computediagnostic!(layer::Layer, ps::CoupledProcesses, state) = _invoke_sequential(computediagnostic!, layer, ps, state)

CryoGrid.computefluxes!(layer::Layer, ps::CoupledProcesses, state) = _invoke_sequential(computefluxes!, layer, ps, state)

CryoGrid.resetfluxes!(layer::Layer, ps::CoupledProcesses, state) = _invoke_sequential(resetfluxes!, layer, ps, state)

function CryoGrid.caninteract(l1::Layer, p1::CoupledProcesses, l2::Layer, p2::Process, state1, state2)
    return any(
        Utils.fastmap(p1.processes) do p_i
            caninteract(l1, p_i, l2, p2, state1, state2)
        end
    )
end

function CryoGrid.caninteract(l1::Layer, p1::Process, l2::Layer, p2::CoupledProcesses, state1, state2)
    return any(
        Utils.fastmap(p2.processes) do p_i
            caninteract(l1, p1, l2, p_i, state1, state2)
        end
    )
end

function CryoGrid.caninteract(l1::Layer, p1::CoupledProcesses, l2::Layer, p2::CoupledProcesses, state1, state2)
    return any(
        Utils.fastmap(p1.processes) do p_i
            nested_caninteract = Utils.fastmap(p2.processes) do p_j
                CryoGrid.caninteract(l1, p_i, l2, p_j, state1, state2)
            end
            return any(nested_caninteract)
        end
    )
end

function CryoGrid.interact!(l1::Layer, p1::CoupledProcesses, l2::Layer, p2::Process, state1, state2)
    Utils.fastiterate(p1.processes) do p_i
        CryoGrid.interactmaybe!(l1, p_i, l2, p2, state1, state2)
    end
end

function CryoGrid.interact!(l1::Layer, p1::Process, l2::Layer, p2::CoupledProcesses, state1, state2)
    Utils.fastiterate(p2.processes) do p_i
        CryoGrid.interactmaybe!(l1, p1, l2, p_i, state1, state2)
    end
end

function CryoGrid.interact!(l1::Layer, p1::CoupledProcesses, l2::Layer, p2::CoupledProcesses, state1, state2)
    Utils.fastiterate(p1.processes) do p_i
        Utils.fastiterate(p2.processes) do p_j
            CryoGrid.interactmaybe!(l1, p_i, l2, p_j, state1, state2)
        end
    end
end

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
