struct TestGroundLayer <: SubSurface end
struct TestGroundProcess <: SubSurfaceProcess end
struct TestBoundary <: BoundaryProcess{TestGroundProcess} end
struct DummyInitializer{varname} <: CryoGrid.Numerics.VarInitializer{varname} end
