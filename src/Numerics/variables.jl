# Variable dimensions
abstract type VarDim{S} end
struct OnGrid{S,F} <: VarDim{S}
    f::F # G -> G' where G is the edge grid
    OnGrid(::Type{S}, f::Function=identity) where {S<:GridSpec} = new{S,typeof(f)}(f)
end
struct Shape{S} <: VarDim{S} Shape(dims::Int...) = new{dims}() end
const Scalar = Shape()
ConstructionBase.constructorof(::Type{OnGrid{S}}) where {S} = f -> OnGrid(S, f)
ConstructionBase.constructorof(::Type{Shape{S}}) where {S} = f -> Shpae(S...)
dimlength(::Shape{dims}, grid::Grid) where dims = prod(dims)
dimlength(d::OnGrid{Cells}, grid::Grid) = d.f(length(cells(grid)))
dimlength(d::OnGrid{Edges}, grid::Grid) = d.f(length(edges(grid)))

abstract type Var{name,S<:VarDim,T,units} end
"""
    Prognostic{name,S,T,units} <: Var{name,S,T,units}

Defines a prognostic (time-integrated) state variable.
"""
struct Prognostic{name,S,T,units} <: Var{name,S,T,units}
    dim::S
    Prognostic(name::Symbol, dims::Union{<:Shape,OnGrid{Cells,typeof(identity)}}, units=NoUnits, ::Type{T}=Float64) where {T} = new{name,typeof(dims),T,units}(dims)
    Prognostic(::Symbol, dims::OnGrid, args...) where {T} = error("Off-cell prognostic/algebraic spatial variables are not currently supported.")
end
"""
    Algebraic{name,S,T,units} <: Var{name,S,T,units}

Defines an algebraic (implicit) state variable.
"""
struct Algebraic{name,S,T,units} <: Var{name,S,T,units}
    dim::S
    # maybe a mass matrix init function?
    Algebraic(name::Symbol, dims::Union{<:Shape,OnGrid{Cells,typeof(identity)}}, units=NoUnits, ::Type{T}=Float64) where {T} = new{name,typeof(dims),T,units}(dims)
    Algebraic(::Symbol, dims::OnGrid, args...) where {T} = error("Off-cell prognostic/algebraic spatial variables are not currently supported.")
end
"""
    Delta{dname,name,S,T,units} <: Var{dname,S,T,units}

Defines a "delta" term `du` for variable `u`, which is the time-derivative or flux for prognostic variables and
the residual for algebraic variables.
"""
struct Delta{dname,name,S,T,units} <: Var{dname,S,T,units}
    dim::S
    Delta(::Symbol, dims::OnGrid, args...) where {T} = error("Off-cell prognostic/algebraic spatial variables are not currently supported.")
    Delta(dname::Symbol, name::Symbol, dims::Union{<:Shape,OnGrid{Cells,typeof(identity)}}, units=NoUnits, ::Type{T}=Float64) where {T} = new{dname,name,typeof(dims),T,units}(dims)
    Delta(var::Prognostic{name,S,T,units}) where {name,S,T,units} = let dims=vardims(var); new{Symbol(:d,name),name,typeof(dims),T,upreferred(units)/u"s"}(dims) end
    Delta(var::Algebraic{name,S,T,units}) where {name,S,T,units} = let dims=vardims(var); new{Symbol(:d,name),name,typeof(dims),T,units}(dims) end
end
"""
    Diagnostic{name,S,T,units} <: Var{name,S,T,units}

Defines a diagnostic variable which is allocated and cached per timestep but not integrated over time.
"""
struct Diagnostic{name,S,T,units} <: Var{name,S,T,units}
    dim::S
    Diagnostic(name::Symbol, dims::VarDim, units=NoUnits, ::Type{T}=Float64) where {T} = new{name,typeof(dims),T,units}(dims)
end
# constructors for Flatten.reconstruct
ConstructionBase.constructorof(::Type{Prognostic{name,S,T,units}}) where {name,S,T,units} = s -> Prognostic(name, s, T, units)
ConstructionBase.constructorof(::Type{Diagnostic{name,S,T,units}}) where {name,S,T,units} = s -> Diagnostic(name, s, T, units)
ConstructionBase.constructorof(::Type{Algebraic{name,S,T,units}}) where {name,S,T,units} = s -> Algebraic(name, s, T, units)
ConstructionBase.constructorof(::Type{Delta{dname,name,S,T,units}}) where {dname,name,S,T,units} = s -> Delta(dname, name, s, T, units)
==(var1::Var{N1,S1,T1,u1},var2::Var{N2,S2,T2,u2}) where {N1,N2,S1,S2,T1,T2,u1,u2} = (N1==N2) && (S1==S2) && (T1==T2) && (u1 == u2)
varname(::Var{name}) where {name} = name
varname(::Type{<:Var{name}}) where {name} = name
vartype(::Var{name,S,T}) where {name,S,T} = T
vartype(::Type{<:Var{name,S,T}}) where {name,S,T} = T
varunits(::Var{name,S,T,units}) where {name,S,T,units} = units
varunits(::Type{<:Var{name,S,T,units}}) where {name,S,T,units} = units
vardims(var::Var) = var.dim
isprognostic(::T) where {T<:Var} = T <: Prognostic
isalgebraic(::T) where {T<:Var} = T <: Algebraic
isflux(::T) where {T<:Var} = T <: Delta
isdiagnostic(::T) where {T<:Var} = T <: Diagnostic
isongrid(::Var{name,S}) where {name,S} = S <: OnGrid
