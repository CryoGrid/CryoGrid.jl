struct Profile{N,D<:DistQuantity,T}
    values::NTuple{N,Pair{D,T}}
    Profile(values::NTuple{N,Pair{D,T}}) where {N,D,T} = new{N,D,T}(values)
    Profile(values::Pair{D,T}...) where {D,T} = new{length(values),D,T}(values)
end
Flatten.flattenable(::Type{<:Profile}, ::Type{Val{:values}}) = false
Base.length(::Profile{N}) where N = N
Base.iterate(profile::Profile) = iterate(profile.values)
Base.iterate(profile::Profile, state) = iterate(profile.values, state)
StructTypes.StructType(::Type{<:Profile}) = StructTypes.CustomStruct()
StructTypes.lower(profile::Profile) = [(depth=StructTypes.lower(row[1]), value=row[2]) for row in profile.values]
StructTypes.lowertype(::Type{<:Profile{N,D,T}}) where {N,D,T} = Vector{NamedTuple{(:depth,:value),Tuple{StructTypes.lowertype(D),Vector{StructTypes.lowertype(T)}}}}
function StructTypes.construct(::Type{<:Profile{N,D,T}}, values::Vector) where {N,D,T}
    depths = [StructTypes.construct(D, row["depth"]) for row in values]
    values = [StructTypes.construct.(T, row["value"]...) for row in values]
    sortinds = sortperm(depths)
    return Profile(collect(map((d,v) -> d => v, depths[sortinds], values[sortinds])))
end

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

"""
    interpolateprofile!(profilearr::DimArray, state; interp=Linear(), extrap=Flat())

Interpolates the given profile to the corresponding variable grids. Assumes state to be indexable via the corresponding
variable symbol and that the parameter names in state and profile match.
"""
function interpolateprofile!(profilearr::DimArray, state; interp=Linear(), extrap=Flat())
    let (depths,names) = dims(profilearr),
        z = ustrip.(depths);
        for p in names
            f = extrapolate(interpolate((z,), ustrip.(profilearr[:,At(p)]), Gridded(interp)), extrap)
            state[p] .= f.(state.grids[p])   # assume length(grid) == length(state.p)
        end
    end
end