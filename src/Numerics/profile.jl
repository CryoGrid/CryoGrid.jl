"""
    Profile{N,idxTypes,valTypes}

Represents a "profile" of values indexed typically over depth, i.e:

x₁: value 1
x₂: value 2
...

where xᵢ are the indices.
"""
struct Profile{N,idxTypes,valTypes}
    indices::idxTypes
    values::valTypes
    Profile(::Tuple{}) = new{0,Tuple{},Tuple{}}((),())
    Profile(indices::NTuple{N,Any}, values::NTuple{N,Any}) where {N} = new{N,typeof(indices),typeof(values)}(indices, values)
    Profile(pairs::Tuple{Vararg{Pair}}) = Profile(map(first, pairs), map(last, pairs))
    Profile(pairs::Pair...) = Profile(pairs)
end
function Base.show(io::IO, mime::MIME"text/plain", profile::Profile)
    for (d,v) in profile
        print("$d ")
        show(io, mime, v)
        println(io)
    end
end
Base.length(::Profile{N}) where N = N
Base.keys(profile::Profile) = profile.indices
Base.values(profile::Profile) = profile.values
Base.iterate(profile::Profile) = iterate(map(Pair, profile.indices, profile.values))
Base.iterate(profile::Profile, state) = iterate(map(Pair, profile.indices, profile.values), state)
Base.map(f, profile::Profile) = Profile(map((i,v) -> i => f(v), profile.indices, profile.values))
Base.getindex(profile::Profile, itrv::Interval) = Profile(Tuple(i => v for (i,v) in profile if i ∈ itrv))
Base.getindex(profile::Profile, i::Int) = profile.indices[i] => profile.values[i]
Base.getindex(profile::Profile, i) = Profile(profile.indices[i], profile.values[i])
Base.lastindex(profile::Profile) = lastindex(profile.indices)
StructTypes.StructType(::Type{<:Profile}) = StructTypes.UnorderedStruct()