import CryoGrid

struct TestGroundLayer{TProc} <: CryoGrid.SubSurface
    proc::TProc
end
CryoGrid.processes(layer::TestGroundLayer) = layer.proc
struct TestGroundProcess <: CryoGrid.SubSurfaceProcess end
struct TestBoundary <: CryoGrid.BoundaryProcess{TestGroundProcess} end
struct DummyInitializer{varname} <: CryoGrid.Strat.VarInitializer{varname} end
