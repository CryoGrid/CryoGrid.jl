# Boundary conditions
struct Snowfall{Tsn<:Forcing{u"m/s"}} <: BoundaryProcess{SnowMassBalance}
    snowfall::Tsn
end
CryoGrid.BoundaryCondition(::Snowfall) = CryoGrid.Neumann()
@inline boundaryvalue(bc::Snowfall, ::Top, ::SnowMassBalance, ::Snowpack, s1, s2) = bc.snowfall(s1.t)
