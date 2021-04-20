"""
      Forcing{T,N}
"""
abstract type Forcing{T,N} end
"""
      TimeSeriesForcing{T,A}
"""
struct TimeSeriesForcing{T,A,I} <: Forcing{T,1}
      tarray::TimeArray{T,1,DateTime,A}
      interp::I
      function TimeSeriesForcing(values::A, timestamps::AbstractArray{DateTime,1}, name::Symbol; interp=Linear()) where
            {T,A<:AbstractArray{T,1}}
            tarray = TimeArray(timestamps,values,[name])
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

export TimeSeriesForcing, Forcing

"""
Get interpolated forcing value at t seconds from t0.
"""
(forcing::TimeSeriesForcing)(t::Real) = forcing.interp(t) # extract interpolation and evaluate
(forcing::TimeSeriesForcing)(t::DateTime) = forcing(ustrip(u"s", float(Dates.datetime2epochms(t))u"ms"))

Base.getindex(f::TimeSeriesForcing, i) = forcing.tarray[i]
function Base.getindex(f::TimeSeriesForcing{T,A,I}, range::StepRange{DateTime,TStep}) where {T,A,I,TStep}
      order(::Interpolations.GriddedInterpolation{T1,N,T2,Gridded{Torder}}) where {T1,N,T2,Torder} = Torder()
      subseries = f.tarray[range]
      TimeSeriesForcing(values(subseries), timestamp(subseries), colnames(f.tarray)[1]; interp=order(f.interp))
end
