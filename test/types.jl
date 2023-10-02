import CryoGrid
import CryoGrid.Hydrology

struct DummyInitializer{varname} <: CryoGrid.VarInitializer{varname} end

struct TestGroundProcess <: CryoGrid.SubSurfaceProcess end

struct TestBoundary <: CryoGrid.BoundaryProcess{TestGroundProcess} end

mutable struct TestGroundLayer{TProc} <: CryoGrid.SubSurface
    proc::TProc
    isactive::Bool
    TestGroundLayer(proc::TProc, isactive::Bool=true) where {TProc} = new{TProc}(proc, isactive)
end

CryoGrid.processes(layer::TestGroundLayer) = layer.proc
CryoGrid.isactive(layer::TestGroundLayer) = layer.isactive

Hydrology.maxwater(::TestGroundLayer, ::WaterBalance) = 1.0
Hydrology.minwater(::TestGroundLayer, ::WaterBalance) = 0.01
Hydrology.hydraulicproperties(::TestGroundLayer) = HydraulicProperties()

