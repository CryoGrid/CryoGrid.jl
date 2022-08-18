const Dirichlet = CryoGrid.Dirichlet
const Neumann = CryoGrid.Neumann
# Boundary condition type aliases
const WaterBC = BoundaryProcess{T} where {WaterBalance<:T<:SubSurfaceProcess}
ConstantInfiltration(value::Quantity) = ConstantBC(WaterBalance, Neumann, uconvert(u"m/s", value))
ImpermeableBoundary() = ConstantBC(WaterBalance, Neumann, 0.0u"m/s")
"""
    Rainfall{Train<:Forcing{u"m/s"}} <: BoundaryProcess{WaterBalance}

Basic rainfall boundary condition for `WaterBalance` which simply invokes the given precipitation
forcing at the current time `t`.
"""
struct Rainfall{Train<:Forcing{u"m/s"}} <: BoundaryProcess{WaterBalance}
    rain::Train
end
CryoGrid.BoundaryStyle(::Type{<:Rainfall}) = Neumann()
function CryoGrid.boundaryvalue(bc::Rainfall, ::Top, ::WaterBalance, ::SubSurface, stop, ssub)
    rainfall_rate = bc.rain(stop.t)
    # take the minimum of the current rainfall rate and the hydraulic conductivity at the top of the upper grid cell;
    # note that this assumes rainfall to be in m/s
    return min(rainfall_rate, ssub.kw[1])
end
