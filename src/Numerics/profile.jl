struct ProfileKnot{D<:DistQuantity,T}
    depth::D
    value::T
end
struct Profile{N,D<:DistQuantity,T}
    knots::NTuple{N,ProfileKnot{D,T}}
    Profile(depths::NTuple{N,D}, values::NTuple{N,T}) where {N,D,T} = new{N,D,T}(depths,values)
    Profile(pairs::NTuple{N,Pair{D,T}}) where {N,D,T} = new{N,D,T}(map(first,pairs), map(last,pairs))
    Profile(pairs::Pair{D,T}...) where {D,T} = Profile(pairs)
end
Flatten.flattenable(::Type{<:ProfileKnot}, ::Type{Val{:depth}}) = false
Base.length(::Profile{N}) where N = N
Base.iterate(profile::Profile) = iterate(profile.knots)
Base.iterate(profile::Profile, state) = iterate(profile.knots, state)
Base.getindex(profile::Profile, i::Int) = profile.knots[i]
Base.getindex(profile::Profile, i) = Profile(profile.knots[i])
Base.lastindex(profile::Profile) = lastindex(profile.knots)
StructTypes.StructType(::Type{<:Profile}) = StructTypes.UnorderedStruct()
