"""
    ConstantBC{P,S,T} <: BoundaryProcess{P}

Constant boundary condition (of any type/unit) specified by `value`.
"""
struct ConstantBC{P,S,T} <: BoundaryProcess{P}
    value::T
    ConstantBC(::Type{P}, ::Type{S}, value::T) where {P<:SubSurfaceProcess,S<:BoundaryStyle,T} = new{P,S,T}(value)
end
ConstructionBase.constructorof(::Type{<:ConstantBC{P,S}}) where {P,S} = value -> ConstantBC(P, S, value)
boundaryvalue(bc::ConstantBC,l1,p2,l2,s1,s2) = bc.value

BoundaryStyle(::Type{<:ConstantBC{P,S}}) where {P,S} = S()

"""
    PeriodicBC{P,S,T} <: BoundaryProcess{P}

Periodic boundary condition (of any type/unit) specified by `period`, `amplitude`, and `phaseshift`.
"""
struct PeriodicBC{P,S,T} <: BoundaryProcess{P}
    period::Float"s"
    amplitude::T
    phaseshift::T
    PeriodicBC(::Type{P}, ::Type{S}, period::Q, amplitude::T=one(T), phaseshift::T=one(T)) where
        {P<:SubSurfaceProcess,S<:BoundaryStyle,Q<:TimeQuantity,T} =
        new{P,S,T}(uconvert(u"s",period) |> dustrip, amplitude, phaseshift)
end
ConstructionBase.constructorof(::Type{<:PeriodicBC{P,S}}) where {P,S} = (args...) -> PeriodicBC(P, S, args...)
@inline boundaryvalue(bc::PeriodicBC,l1,p2,l2,s1,s2) = bc.amplitude*sin(π*(1/bc.period)*t + bc.phaseshift)

BoundaryStyle(::Type{<:PeriodicBC{P,S}}) where {P,S} = S()

Base.@kwdef struct Bias{P,B} <: BoundaryProcess{P}
    bias::B = Param(0.0)
    Bias(::Type{P}, bias::B) where {P<:SubSurfaceProcess,B} = new{P,B}(bias)
end
@inline boundaryvalue(bc::Bias,l1,p2,l2,s1,s2) = bc.bias

BoundaryStyle(::Type{<:Bias}) = Dirichlet()

"""
    struct CombinedBoundaryProcess{B1,B2,F,S} <: BoundaryProcess

Represents a composition of two boundary processes, `B1` and `B2`, via an operator `F`.
A typical use case is combining `ConstantBC` with a forcing-driven boundary process to
scale or shift the forcing.
"""
struct CombinedBoundaryProcess{B1,B2,P,F,S} <: BoundaryProcess{P}
    op::F
    bc1::B1
    bc2::B2
    function CombinedBoundaryProcess(op::F, bc1::B1, bc2::B2) where {F,P,B1<:BoundaryProcess{P},B2<:BoundaryProcess{P}}
        @assert BoundaryStyle(bc1) == BoundaryStyle(bc2) "boundary condition styles (e.g. Dirichlet vs Neumann) must match"
        new{B1,B2,P,F,typeof(BoundaryStyle(bc1))}(op,bc1,bc2)
    end
end
@inline boundaryvalue(cbc::CombinedBoundaryProcess{B1,B2},l1,p2,l2,s1,s2) where {B1<:BoundaryProcess,B2<:BoundaryProcess} = cbc.op(boundaryvalue(cbc.bc1,l1,l2,p2,s1,s2), boundaryvalue(cbc.bc2,l1,l2,p2,s1,s2))
@inline boundaryvalue(cbc::CombinedBoundaryProcess{B1,B2},l1,p2,l2,s1,s2) where {B1,B2<:BoundaryProcess} = cbc.op(cbc.bc1, boundaryvalue(cbc.bc2,l1,l2,p2,s1,s2))
@inline boundaryvalue(cbc::CombinedBoundaryProcess{B1,B2},l1,p2,l2,s1,s2) where {B1<:BoundaryProcess,B2} = cbc.op(boundaryvalue(cbc.bc1,l1,l2,p2,s1,s2), cbc.bc2)
variables(top::Top, cbc::CombinedBoundaryProcess) = tuplejoin(variables(top, cbc.bc1), variables(top, cbc.bc2))
BoundaryStyle(::Type{CombinedBoundaryProcess{B1,B2,P,F,S}}) where {F,P,B1,B2,S} = S()
# Overload arithmetic operators on boundary processes.
Base.:+(bc1::BoundaryProcess, bc2::BoundaryProcess) = CombinedBoundaryProcess(+, bc1, bc2)
Base.:-(bc1::BoundaryProcess, bc2::BoundaryProcess) = CombinedBoundaryProcess(-, bc1, bc2)
Base.:*(bc1::BoundaryProcess, bc2::BoundaryProcess) = CombinedBoundaryProcess(*, bc1, bc2)
Base.:/(bc1::BoundaryProcess, bc2::BoundaryProcess) = CombinedBoundaryProcess(/, bc1, bc2)
