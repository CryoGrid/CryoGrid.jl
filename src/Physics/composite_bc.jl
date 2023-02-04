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
        @assert BoundaryStyle(bc1) == BoundaryStyle(bc2) "boundary condition styles (e.g. Dirichlet vs Neumann) must match"
        new{B1,B2,P,F,typeof(BoundaryStyle(bc1))}(op,bc1,bc2)
    end
end
@inline CryoGrid.boundaryvalue(cbc::CompositeBoundaryProcess{B1,B2},l1,p2,l2,s1,s2) where {B1<:BoundaryProcess,B2<:BoundaryProcess} = cbc.op(boundaryvalue(cbc.bc1,l1,p2,l2,s1,s2), boundaryvalue(cbc.bc2,l1,p2,l2,s1,s2))
@inline CryoGrid.boundaryvalue(cbc::CompositeBoundaryProcess{B1,B2},l1,p2,l2,s1,s2) where {B1,B2<:BoundaryProcess} = cbc.op(cbc.bc1, boundaryvalue(cbc.bc2,l1,p2,l2,s1,s2))
@inline CryoGrid.boundaryvalue(cbc::CompositeBoundaryProcess{B1,B2},l1,p2,l2,s1,s2) where {B1<:BoundaryProcess,B2} = cbc.op(boundaryvalue(cbc.bc1,l1,p2,l2,s1,s2), cbc.bc2)
CryoGrid.variables(top::Top, cbc::CompositeBoundaryProcess) = tuplejoin(variables(top, cbc.bc1), variables(top, cbc.bc2))
CryoGrid.BoundaryStyle(::Type{CompositeBoundaryProcess{B1,B2,P,F,S}}) where {F,P,B1,B2,S} = S()
# Overload arithmetic operators on boundary processes.
Base.:+(bc1::BoundaryProcess, bc2::BoundaryProcess) = CompositeBoundaryProcess(+, bc1, bc2)
Base.:-(bc1::BoundaryProcess, bc2::BoundaryProcess) = CompositeBoundaryProcess(-, bc1, bc2)
Base.:*(bc1::BoundaryProcess, bc2::BoundaryProcess) = CompositeBoundaryProcess(*, bc1, bc2)
Base.:/(bc1::BoundaryProcess, bc2::BoundaryProcess) = CompositeBoundaryProcess(/, bc1, bc2)
