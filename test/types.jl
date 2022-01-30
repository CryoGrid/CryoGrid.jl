struct TestGroundLayer <: SubSurface end
struct TestGroundProcess <: SubSurfaceProcess end
struct TestBoundary <: BoundaryProcess end
struct DummyInitializer{varname} <: CryoGrid.Numerics.VarInitializer{varname} end
