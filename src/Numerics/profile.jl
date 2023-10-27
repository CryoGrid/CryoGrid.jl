struct ProfileKnot{D,V}
    depth::D
    value::V
end
Base.iterate(knot::ProfileKnot) = iterate((knot.depth, knot.value))
Base.iterate(knot::ProfileKnot, state) = iterate((knot.depth, knot.value), state)
Base.show(io::IO, knot::ProfileKnot) = print(io, "$(knot.depth): $(knot.value)")
struct Profile{N,V,D}
    knots::NTuple{N,ProfileKnot{D,V}}
    Profile(::Tuple{}) = new{0,Nothing,Nothing}(())
    Profile(knots::Tuple{Vararg{ProfileKnot{D,V},N}}) where {N,D,V} = new{N,V,D}(knots)
    Profile(pairs::Tuple{Vararg{Pair}}) = Profile(map(Base.splat(ProfileKnot), pairs))
    Profile(pairs::Pair...) = Profile(pairs)
end
function Base.show(io::IO, mime::MIME"text/plain", profile::Profile)
    for knot in profile
        show(io, mime, knot)
        println(io)
    end
end
Base.length(::Profile{N}) where N = N
Base.iterate(profile::Profile) = iterate(profile.knots)
Base.iterate(profile::Profile, state) = iterate(profile.knots, state)
Base.map(f, profile::Profile) = Profile(map(knot -> ProfileKnot(knot.depth, f(knot.value)), profile.knots))
Base.getindex(profile::Profile, itrv::Interval) = Profile(Tuple(knot for knot in profile.knots if knot.depth ∈ itrv))
Base.getindex(profile::Profile, i::Int) = profile.knots[i]
Base.getindex(profile::Profile, i) = Profile(profile.knots[i])
Base.lastindex(profile::Profile) = lastindex(profile.knots)
StructTypes.StructType(::Type{<:Profile}) = StructTypes.UnorderedStruct()
