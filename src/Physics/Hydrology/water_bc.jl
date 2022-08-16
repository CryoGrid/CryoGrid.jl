const Dirichlet = CryoGrid.Dirichlet
const Neumann = CryoGrid.Neumann
# Boundary condition type aliases
const WaterBC = BoundaryProcess{T} where {WaterBalance<:T<:SubSurfaceProcess}
ConstantInfiltration(value::Quantity) = ConstantBC(WaterBalance, Neumann, uconvert(u"m/s", value))
ImpermeableBoundary() = ConstantBC(WaterBalance, Neumann, 0.0u"m/s")
"""
    Rainfall{Train<:Forcing} <: BoundaryProcess{WaterBalance}

Basic rainfall boundary condition for `WaterBalance` which simply invokes the given precipitation
forcing at the current time `t`.
"""
struct Rainfall{Train<:Forcing} <: BoundaryProcess{WaterBalance}
    rain::Train
end
CryoGrid.BoundaryStyle(::Type{<:Rainfall}) = Neumann()
@inline CryoGrid.boundaryvalue(bc::Rainfall, l1, ::WaterBalance, l2, s1, s2) = getscalar(s1.jw_ub)
CryoGrid.variables(::Rainfall) = (
    Diagnostic(:jw_ub, Scalar, u"m/s")
)
function CryoGrid.diagnosticstep!(top::Top, bc::Rainfall, state)
    @setscalar state.jw_ub = bc.rain(state.t)
end
