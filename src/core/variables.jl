import Base.==

# Variable dimensions
abstract type VarDim end
struct OnGrid{S} <: VarDim
    f::Function # G -> G' where G is the edge grid
    OnGrid(::Type{S}, f::Function=x->x) where {S<:GridSpec} = new{S}(f)
end
struct Shape{D} <: VarDim Shape(dims::Int...) = new{dims}() end
const Scalar = Shape()

# Variable trait (prognostic, diagnostic, parameter)
abstract type VarStyle end
struct Prognostic <: VarStyle end
struct Diagnostic <: VarStyle end
struct Parameter <: VarStyle end

# Variable type
struct Var{name,T,D,S}
    dim::D
    Var(name::Symbol, typ::Type{T}, dim::Union{<:VarDim,<:GridSpec}=Scalar, style::Type{S}=Diagnostic) where {T,S<:VarStyle} =
        new{name,T,typeof(dim),S}(dim)
end
==(var1::Var{N1,T1,D1},var2::Var{N2,T2,D2}) where {N1,N2,T1,T2,D1,D2} = (N1==N2) && (T1==T2) && (D1==D2)
varname(::Var{name}) where {name} = name
varname(::Type{<:Var{name}}) where {name} = name
vartype(::Var{name,T}) where {name,T} = T
vartype(::Type{<:Var{name,T}}) where {name,T} = T
export Var, varname
export VarDim, OnGrid, Shape, Scalar

Prognostic(name::Symbol, typ::Type{T}, dim::Union{<:VarDim,<:GridSpec}=Scalar) where T = Var(name,typ,dim,Prognostic)
Diagnostic(name::Symbol, typ::Type{T}, dim::Union{<:VarDim,<:GridSpec}=Scalar) where T = Var(name,typ,dim,Diagnostic)
Parameter(name::Symbol, typ::Type{T}=Float64, dim::Union{<:VarDim,<:GridSpec}=Scalar) where T = Var(name,typ,dim,Parameter)

VarStyle(::Type{Var{name,T,D,S}}) where {name,T,D,S} = S()
isprognostic(var::Var{name,T,D,S}) where {name,T,D,S} = S == Prognostic
isdiagnostic(var::Var{name,T,D,S}) where {name,T,D,S} = S == Diagnostic
isparameter(var::Var{name,T,D,S}) where {name,T,D,S} = S == Parameter
export VarStyle, Prognostic, Diagnostic, isprognostic, isdiagnostic, isparameter

struct Vars{ProgVars,NT}
    data::NT
    Vars(progvars::NTuple{N,Symbol}, data::NamedTuple) where {N} = new{progvars,typeof(data)}(data)
end
# Passthrough to named tuple
Base.getproperty(vars::Vars, sym::Symbol) = getproperty(getfield(vars, :data), sym)
