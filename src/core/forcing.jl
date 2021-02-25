"""
      Forcing{T,N}
"""
abstract type Forcing{T,N} end
"""
      TimeSeriesForcing{T,A}
"""
struct TimeSeriesForcing{T,A} <: Forcing{T,1}
      t0::DateTime
      tarray::TimeArray{T,1,DateTime,A}
      TimeSeriesForcing(values::A, timestamps::AbstractArray{DateTime,1}, name::Symbol; interp=Linear()) where
            {T,A<:AbstractArray{T,1}} = begin
            t0 = timestamps[1]
            ts = Dates.value.(timestamps .- t0)./1000.0
            interp = interpolate((ts,), values, Gridded(interp))
            new{T,typeof(interp)}(t0,TimeArray(timestamps, interp, [name]))
      end
end

export TimeSeriesForcing, Forcing

"""
Get forcing at time t.
"""
(forcing::TimeSeriesForcing)(t) = values(forcing.tarray)(t) # extract interpolation and evaluate
