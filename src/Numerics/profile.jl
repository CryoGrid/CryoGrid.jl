struct Profile{N,D<:DistQuantity,T}
    depths::NTuple{N,D}
    values::NTuple{N,T}
    Profile(depths::NTuple{N,D}, values::NTuple{N,T}) where {N,D,T} = new{N,D,T}(depths,values)
    Profile(pairs::NTuple{N,Pair{D,T}}) where {N,D,T} = new{N,D,T}(map(first,pairs), map(last,pairs))
    Profile(pairs::Pair{D,T}...) where {D,T} = Profile(pairs)
end
Flatten.flattenable(::Type{<:Profile}, ::Type{Val{:depths}}) = false
Base.length(::Profile{N}) where N = N
Base.iterate(profile::Profile) = iterate(zip(profile.depths,profile.values))
Base.iterate(profile::Profile, state) = iterate(zip(profile.depths,profile.values),state)
Base.getindex(profile::Profile, i) = (depth=profile.depths[i], value=profile.values[i])
StructTypes.StructType(::Type{<:Profile}) = StructTypes.UnorderedStruct()
"""
    profile2array(profile::Profile{N,D,T};names) where {N,D,T}

Constructs a DimArray from the given Profile, i.e. pairs Q => (x1,...,xn) where x1...xn are the values defined at Q.
Column names for the resulting DimArray can be set via the names parameter which accepts an NTuple of symbols,
where N must match the number of parameters given (i.e. n).
"""
function profile2array(profile::Profile{N,D,T};names::Union{Nothing,NTuple{M,Symbol}}=nothing) where {M,N,D,T}
    depths, vals = zip(profile.values...)
    params = hcat(collect.(vals)...)'
    names = isnothing(names) ? [Symbol(:x,:($i)) for i in 1:N] : collect(names)
    DimArray(params, (Z(collect(depths)), Dim{:var}(names)))
end
array2profile(arr::AbstractDimArray) = Profile(collect(dims(profile,Z)), mapslices(collect, arr; dims=2))
