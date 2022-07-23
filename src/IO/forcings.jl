"""
      Forcing{T,N}

Abstract type representing a generic external boundary condition (i.e. "forcing").
"""
abstract type Forcing{T,N} end
@inline @propagate_inbounds (forcing::Forcing)(x::Number) = error("$(typeof(forcing)) not implemented")
@inline @propagate_inbounds (forcing::Forcing)(t::DateTime) = forcing(ustrip(u"s", float(Dates.datetime2epochms(t))u"ms"))
# disable flattening for all fields of forcing types by default
Flatten.flattenable(::Type{<:Forcing}, ::Type) = false
InputOutput.parameterize(f::Forcing; fields...) = f

"""
      TimeSeriesForcing{T,A,I}

Forcing provided by a discrete time series of data.
"""
struct TimeSeriesForcing{T,A,I} <: Forcing{T,1}
      tarray::TimeArray{T,1,DateTime,A}
      interp::I
      TimeSeriesForcing(tarray::TimeArray{T,N,D,A}, interp) where {T,N,D,A} = new{T,A,typeof(interp)}(tarray,interp)
      function TimeSeriesForcing(values::A, timestamps::AbstractArray{DateTime,1}, name::Symbol; interp=Linear()) where
            {T,A<:AbstractArray{T,1}}
            tarray = TimeArray(Vector(timestamps),Vector(values),[name])
            TimeSeriesForcing(tarray; interp=interp)
      end
      function TimeSeriesForcing(tarray::TimeArray{T,N,D,A}; interp=Linear()) where {T,N,D,A}
            # converts timestamps to seconds since epoch
            ts = ustrip.(u"s", float.(Dates.datetime2epochms.(timestamp(tarray)))u"ms")
            interpolated = interpolate((ts,), values(tarray), Gridded(interp))
            new{T,A,typeof(interpolated)}(tarray,interpolated)
      end
end

Base.show(io::IO, forcing::TimeSeriesForcing{T}) where T = print(io, "TimeSeriesForcing{$T}($(typeof(forcing.tarray))}")

"""
Get interpolated forcing value at t seconds from t0.
"""
@inline @propagate_inbounds (forcing::TimeSeriesForcing)(t::Number) = forcing.interp(t) # extract interpolation and evaluate

@inline @propagate_inbounds Base.getindex(forcing::TimeSeriesForcing, i) = forcing.tarray[i]
function Base.getindex(f::TimeSeriesForcing{T,A,I}, range::StepRange{DateTime,TStep}) where {T,A,I,TStep}
      order(::Interpolations.GriddedInterpolation{T1,N,T2,Gridded{Torder}}) where {T1,N,T2,Torder} = Torder()
      subseries = f.tarray[range]
      TimeSeriesForcing(values(subseries), timestamp(subseries), colnames(f.tarray)[1]; interp=order(f.interp))
end
