"""
    struct CompositeBoundaryProcess{F,B1,B2} <: BoundaryProcess

Represents a composition of two boundary processes, `B1` and `B2`, via an operator `F`.
A typical use case is combining `Constant` with a forcing-driven boundary process to
scale or shift the forcing.
"""
struct CompositeBoundaryProcess{F,B1,B2,S} <: BoundaryProcess
    op::F
    bc1::B1
    bc2::B2
    function CompositeBoundaryProcess(op::F, bc1::B1, bc2::B2) where {F,B1<:BoundaryProcess,B2<:BoundaryProcess}
        @assert BoundaryStyle(bc1) == BoundaryStyle(bc2) "boundary condition styles (e.g. Dirichlet vs Neumann) must match"
        new{F,B1,B2,typeof(BoundaryStyle(bc1))}(op,bc1,bc2)
    end
end

@inline (cbc::CompositeBoundaryProcess)(l1,l2,p2,s1,s2) = cbc.op(cbc.bc1(l1,l2,p2,s1,s2), cbc.bc2(l1,l2,p2,s1,s2))

variables(top::Top, cbc::CompositeBoundaryProcess) = tuplejoin(variables(top, cbc.bc1), variables(top, cbc.bc2))
BoundaryStyle(::Type{CompositeBoundaryProcess{F,B1,B2,S}}) where {F,B1,B2,S} = S()

# Overload arithmetic operators on boundary processes.
Base.:+(bc1::BoundaryProcess, bc2::BoundaryProcess) = CompositeBoundaryProcess(+, bc1, bc2)
Base.:-(bc1::BoundaryProcess, bc2::BoundaryProcess) = CompositeBoundaryProcess(-, bc1, bc2)
Base.:*(bc1::BoundaryProcess, bc2::BoundaryProcess) = CompositeBoundaryProcess(*, bc1, bc2)
Base.:/(bc1::BoundaryProcess, bc2::BoundaryProcess) = CompositeBoundaryProcess(/, bc1, bc2)
