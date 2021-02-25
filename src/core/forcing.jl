"""
      Forcing{T,N}
"""
abstract type Forcing{T,N} end
"""
      TimeSeriesForcing{T,A}
"""
struct TimeSeriesForcing{T,A,I} <: Forcing{T,1}
      t0::DateTime
      tarray::TimeArray{T,1,DateTime,A}
      interp::I
      TimeSeriesForcing(values::A, timestamps::AbstractArray{DateTime,1}, name::Symbol; interp=Linear()) where
            {T,A<:AbstractArray{T,1}} = begin
            t0 = timestamps[1]
            ts = Dates.value.(timestamps .- t0)./1000.0
            tarray = TimeArray(timestamps,values,[name])
            interpolated = interpolate((ts,), values, Gridded(interp))
            new{T,A,typeof(interpolated)}(t0,tarray,interpolated)
      end
end

export TimeSeriesForcing, Forcing

"""
Get interpolated forcing value at t seconds from t0.
"""
(forcing::TimeSeriesForcing)(t) = forcing.interp(t) # extract interpolation and evaluate

Base.getindex(f::TimeSeriesForcing, i) = forcing.tarray[i]
function Base.getindex(f::TimeSeriesForcing{T,A,I}, range::StepRange{DateTime,TStep}) where {T,A,I,TStep}
      order(::Interpolations.GriddedInterpolation{T1,N,T2,Gridded{Torder}}) where {T1,N,T2,Torder} = Torder()
      subseries = f.tarray[range]
      TimeSeriesForcing(values(subseries), timestamp(subseries), colnames(f.tarray)[1]; interp=order(f.interp))
end
