"""
    CompositeBoundaryProcess{B1,B2,F,S} <: BoundaryProcess

Represents a composition of two boundary processes, `B1` and `B2`, via an operator `F`.
A typical use case is combining `ConstantBC` with a forcing-driven boundary process to
scale or shift the forcing.
"""
struct CompositeBoundaryProcess{B1,B2,P,F,S} <: BoundaryProcess{P}
    op::F
    bc1::B1
    bc2::B2
    function CompositeBoundaryProcess(op::F, bc1::B1, bc2::B2) where {F,P,B1<:BoundaryProcess{P},B2<:BoundaryProcess{P}}
        @assert BCKind(bc1) == BCKind(bc2) "boundary condition styles (e.g. Dirichlet vs Neumann) must match"
        new{B1,B2,P,F,typeof(BCKind(bc1))}(op,bc1,bc2)
    end
end
CryoGrid.boundaryvalue(cbc::CompositeBoundaryProcess{B1,B2}, state) where {B1<:BoundaryProcess,B2<:BoundaryProcess} = cbc.op(boundaryvalue(cbc.bc1, state), boundaryvalue(cbc.bc2, state))
CryoGrid.variables(top::Top, cbc::CompositeBoundaryProcess) = tuplejoin(variables(top, cbc.bc1), variables(top, cbc.bc2))
CryoGrid.BCKind(::Type{CompositeBoundaryProcess{B1,B2,P,F,S}}) where {F,P,B1,B2,S} = S()
# Overload arithmetic operators on boundary processes.
Base.:+(bc1::BoundaryProcess, bc2::BoundaryProcess) = CompositeBoundaryProcess(+, bc1, bc2)
Base.:-(bc1::BoundaryProcess, bc2::BoundaryProcess) = CompositeBoundaryProcess(-, bc1, bc2)
Base.:*(bc1::BoundaryProcess, bc2::BoundaryProcess) = CompositeBoundaryProcess(*, bc1, bc2)
Base.:/(bc1::BoundaryProcess, bc2::BoundaryProcess) = CompositeBoundaryProcess(/, bc1, bc2)
