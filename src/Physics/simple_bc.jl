"""
    BoundaryEffect

Base type for boundary "effects" which modify boundary conditions based on some
given parameterization.
"""
abstract type BoundaryEffect end

"""
    ConstantBC{P,S,T} <: BoundaryProcess{P}

Constant boundary condition (of any type/unit) specified by `value`.
"""
struct ConstantBC{P,S,T} <: BoundaryProcess{P}
    value::T
    ConstantBC(::Type{P}, ::Type{S}, value::T) where {P<:SubSurfaceProcess,S<:BoundaryStyle,T} = new{P,S,T}(value)
end
ConstructionBase.constructorof(::Type{<:ConstantBC{P,S}}) where {P,S} = value -> ConstantBC(P, S, value)
CryoGrid.boundaryvalue(bc::ConstantBC,l1,p2,l2,s1,s2) = bc.value

CryoGrid.BoundaryStyle(::Type{<:ConstantBC{P,S}}) where {P,S} = S()

"""
    PeriodicBC{P,S,T1,T2,T3,T4} <: BoundaryProcess{P}

Periodic boundary condition (of any type/unit) specified by `period`, `amplitude`, and `phaseshift`.
"""
struct PeriodicBC{P,S,T1,T2,T3,T4} <: BoundaryProcess{P}
    period::T1
    amplitude::T2
    phaseshift::T3
    level::T4
    PeriodicBC(::Type{P}, ::Type{S}, period::T1=1.0, amplitude::T2=1.0, phaseshift::T3=0.0, level::T4=0.0) where
        {P<:SubSurfaceProcess,S<:BoundaryStyle,T1,T2,T3,T4} =
        new{P,S,T1,T2,T3,T4}(period, amplitude, phaseshift, level)
end
ConstructionBase.constructorof(::Type{<:PeriodicBC{P,S}}) where {P,S} = (args...) -> PeriodicBC(P, S, args...)
CryoGrid.boundaryvalue(bc::PeriodicBC,l1,p2,l2,s1,s2) = bc.amplitude*sin(π*(1/bc.period)*s1.t + bc.phaseshift) + bc.level
CryoGrid.BoundaryStyle(::Type{<:PeriodicBC{P,S}}) where {P,S} = S()

# convenience constructors
ConstantValue(::Type{P}, value::T) where {P<:SubSurfaceProcess,T} = ConstantBC(P, Dirichlet, value)
PeriodicValue(::Type{P}, args...) where {P<:SubSurfaceProcess} = PeriodicBC(P, Dirichlet, args...)
ConstantFlux(::Type{P}, value::T) where {P<:SubSurfaceProcess,T} = ConstantBC(P, Neumann, value)
PeriodicFlux(::Type{P}, args...) where {P<:SubSurfaceProcess} = PeriodicBC(P, Neumann, args...)
