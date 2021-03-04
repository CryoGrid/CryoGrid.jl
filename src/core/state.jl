"""
    CryoGridState{T,A,Ax,S} <: DEDataVector{T}

Specialized implementation of DEDataVector that holds non-integrated state variables in a (possibly nested) named
tuple field 'state'. Relevant overrides for AbstractArray and DiffEqBase interface methods are provided to recursively
copy or convert data types of all arrays nested in 'state'.
"""
struct CryoGridState{T,A,Ax,S} <: DEDataVector{T}
    x::A # mandatory field by DEDataArray specification
    ax::Ax # axes for component array
    state::S # diagnostic state
    function CryoGridState(u::A, state::S) where {T,A<:ComponentVector{T},S<:NamedTuple}
        x = getdata(u)
        ax = getaxes(u)
        new{T,typeof(x),typeof(ax),S}(x,ax,state)
    end
    function CryoGridState(x::A, ax::Ax, state::S) where {T,A<:AbstractVector{T},Ax<:Tuple{<:ComponentArrays.Axis},S<:NamedTuple}
        new{T,A,Ax,S}(x,ax,state)
    end
end

withaxes(u::CryoGridState) = ComponentArray(u.x, u.ax)

export CryoGridState, withaxes

# similar CryoGridState
Base.similar(A::CryoGridState{T}) where T = similar(A,T)
Base.similar(A::CryoGridState, ::Type{T}) where T = CryoGridState(similar(A.x,T),A.ax,copystate(A.state,T))
Base.similar(A::CryoGridState, ::Type{T}, dims::NTuple{N,Int}) where {T,N} = similar(A.x,T,dims)
# copy state named tuples (nested tuples of arrays);
# recurses through each named tuple calling simlar on all arrays.
copystate(s::NamedTuple, ::Type{T}) where T = map(x -> copystate(x,T), s)
# base case 1: similar array type
copystate(x::AbstractArray, ::Type{T}) where T = copyto!(similar(x,T),x)
# base case 2: grid type
copystate(x::Grid, ::Type{T}) where T = x
# base case 3: any other type, just copy (assume immutable)
copystate(x, ::Type{T}) where T = x
# DiffEqBase copy_fields
DiffEqBase.copy_fields!(arr::CryoGridState{TDest}, template::CryoGridState{TSrc}) where {TDest,TSrc} =
    error("mutating copy_fields! is not supported for immutable CryoGridState")
DiffEqBase.copy_fields(arr::AbstractArray{TDest}, template::CryoGridState{TSrc}) where {TDest,TSrc} =
    CryoGridState(arr,template.ax,copystate(template.state,TDest))
# recursive array copy
RecursiveArrayTools.recursivecopy!(dest::T, src::T) where {T<:Tuple} = src
RecursiveArrayTools.recursivecopy!(dest::T, src::T) where {T<:Grid} = src
RecursiveArrayTools.recursivecopy!(dest::T, src::T) where {T<:NamedTuple} = map(x -> recursivecopy!(x[2],x[1]), zip(src,dest))
function RecursiveArrayTools.recursivecopy!(dest::T, src::T) where {T<:CryoGridState}
    @inbounds copyto!(dest,src)
    @inbounds recursivecopy!(dest.state, src.state)
    return dest
end
