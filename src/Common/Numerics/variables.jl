# Variable dimensions
abstract type VarDim end
struct OnGrid{S} <: VarDim
    f::Function # G -> G' where G is the edge grid
    OnGrid(::Type{S}, f::Function=identity) where {S<:GridSpec} = new{S}(f)
end
struct Shape{S} <: VarDim Shape(dims::Int...) = new{dims}() end
const Scalar = Shape()

abstract type Var{name,T,S<:Union{<:VarDim,<:GridSpec}} end
struct Prognostic{name,T,S} <: Var{name,T,S}
    dim::S
    Prognostic(name::Symbol, ::Type{T}, dims::Union{<:VarDim,<:GridSpec}) where {T} = new{name, T, typeof(dims)}(dims)
end
struct Algebraic{name,T,S} <: Var{name,T,S}
    dim::S
    # maybe a mass matrix init function?
    Algebraic(name::Symbol, ::Type{T}, dims::Union{<:VarDim,<:GridSpec}) where {T} = new{name, T, typeof(dims)}(dims)
end
struct Diagnostic{name,T,S} <: Var{name,T,S}
    dim::S
    Diagnostic(name::Symbol, ::Type{T}, dims::Union{<:VarDim,<:GridSpec}) where {T} = new{name, T, typeof(dims)}(dims)
end
struct Parameter{name,T,S,D} <: Var{name,T,S}
    dim::S
    default_value::T
    Parameter(name::Symbol, default_value::T, domain::Interval{L,R,I}=-Inf..Inf) where {L,R,I<:Real,T<:AbstractArray} =
        new{name,T,typeof(Shape(size(default_value)...)),convert(Interval{L,R,Float64}, domain)}(Shape(size(default_value)...), default_value)
    Parameter(name::Symbol, default_value::T, domain::Interval=-Inf..Inf) where {T<:Real} = Parameter(name, [default_value], domain)
end
ConstructionBase.constructorof(::Type{Prognostic{name,T,S}}) where {name,T,S} = s -> Prognostic(name, T, s)
ConstructionBase.constructorof(::Type{Diagnostic{name,T,S}}) where {name,T,S} = s -> Diagnostic(name, T, s)
ConstructionBase.constructorof(::Type{Algebraic{name,T,S}}) where {name,T,S} = s -> Algebraic(name, T, s)
==(var1::Var{N1,T1,D1},var2::Var{N2,T2,D2}) where {N1,N2,T1,T2,D1,D2} = (N1==N2) && (T1==T2) && (D1==D2)
varname(::Var{name}) where {name} = name
varname(::Type{<:Var{name}}) where {name} = name
vartype(::Var{name,T}) where {name,T} = T
vartype(::Type{<:Var{name,T}}) where {name,T} = T
vardims(var::Var{name,T,S}) where {name,T,S} = var.dim
isprognostic(::T) where {T<:Var} = T <: Prognostic
isalgebraic(::T) where {T<:Var} = T <: Algebraic
isdiagnostic(::T) where {T<:Var} = T <: Diagnostic
isparameter(::T) where {T<:Var} = T <: Parameter
domain(::Parameter{name,T,S,D}) where {name,T,S,D} = D

# parameter constraints
function checkdomain(domain::Interval, f::Function, z)
    @assert all(z .∈ 0..Inf) "value $z for parameter $name is outside of the given domain: $domain"
    f.(z)
end
constrain(::Parameter{name,T,S,-Inf..Inf}, x) where {name,T,S} = x
unconstrain(::Parameter{name,T,S,-Inf..Inf}, z) where {name,T,S} = z
constrain(::Parameter{name,T,S,0.0..1.0}, x) where {name,T,S} = logistic.(x)
unconstrain(::Parameter{name,T,S,0.0..1.0}, z) where {name,T,S} = checkdomain(0.0..1.0, logit, z)
constrain(::Parameter{name,T,S,0..Inf}, x) where {name,T,S} = softplus.(x)
unconstrain(::Parameter{name,T,S,0..Inf}, z) where {name,T,S} = checkdomain(0..Inf, softplusinv, z)
constrain(::Parameter{name,T,S,1..Inf}, x) where {name,T,S} = (plusone ∘ softplus).(x)
unconstrain(::Parameter{name,T,S,1..Inf}, z) where {name,T,S} = checkdomain(1..Inf, softplusinv ∘ minusone, z)
