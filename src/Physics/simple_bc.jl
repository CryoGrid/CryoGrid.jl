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
    ConstantBC(::Type{P}, ::Type{S}, value::T) where {P<:SubSurfaceProcess,S<:BCKind,T} = new{P,S,T}(value)
end
ConstructionBase.constructorof(::Type{<:ConstantBC{P,S}}) where {P,S} = value -> ConstantBC(P, S, value)
CryoGrid.boundaryvalue(bc::ConstantBC, state) = bc.value
CryoGrid.BCKind(::Type{<:ConstantBC{P,S}}) where {P,S} = S()

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
        {P<:SubSurfaceProcess,S<:BCKind,T1,T2,T3,T4} =
        new{P,S,T1,T2,T3,T4}(period, amplitude, phaseshift, level)
end
ConstructionBase.constructorof(::Type{<:PeriodicBC{P,S}}) where {P,S} = (args...) -> PeriodicBC(P, S, args...)
CryoGrid.boundaryvalue(bc::PeriodicBC, state) = bc.amplitude*sin(2π*(state.t/bc.period) + bc.phaseshift) + bc.level
CryoGrid.BCKind(::Type{<:PeriodicBC{P,S}}) where {P,S} = S()

# convenience constructors
ConstantValue(::Type{P}, value::T) where {P<:SubSurfaceProcess,T} = ConstantBC(P, Dirichlet, value)
PeriodicValue(::Type{P}, args...) where {P<:SubSurfaceProcess} = PeriodicBC(P, Dirichlet, args...)
ConstantFlux(::Type{P}, value::T) where {P<:SubSurfaceProcess,T} = ConstantBC(P, Neumann, value)
PeriodicFlux(::Type{P}, args...) where {P<:SubSurfaceProcess} = PeriodicBC(P, Neumann, args...)

"""
    FunctionBC{P,S,F} <: BoundaryProcess{P}

Generic boundary condition type that invokes a user supplied function `func(state)` where `state` is the
layer state associated with the boundary layer.
"""
struct FunctionBC{P,S,F} <: BoundaryProcess{P}
    func::F
    FunctionBC(::Type{P}, ::Type{S}, func::F) where {P<:SubSurfaceProcess,S<:BCKind,F} = new{P,S,F}(func)
end
ConstructionBase.constructorof(::Type{<:FunctionBC{P,S}}) where {P,S} = f -> FunctionBC(P, S, f)
CryoGrid.boundaryvalue(bc::FunctionBC, state) = bc.func(state)
CryoGrid.BCKind(::Type{<:FunctionBC{P,S}}) where {P,S} = S()

