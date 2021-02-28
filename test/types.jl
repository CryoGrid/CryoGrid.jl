struct TestGroundLayer <: SubSurface end
struct TestGroundProcess <: SubSurfaceProcess end
struct TestBoundary <: BoundaryProcess{TestGroundProcess} end

const DummySoilProfile = SoilProfile(
	0.0u"m" => (0.5,0.0,0.5,0.0,0.0),
	1000.0u"m" => (0.5,0.0,0.5,0.0,0.0),
)

const DummyTempProfile = TempProfile(
	0.0u"m" => 0.0u"°C",
	1000.0u"m" => 0.0u"°C",
)
