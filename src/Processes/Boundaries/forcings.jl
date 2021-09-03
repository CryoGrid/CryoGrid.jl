ForcingData = (
      Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044 = Resource("Samoylov_ERA_obs_fitted_1979_2014_spinup_extended_2044", "json", "https://nextcloud.awi.de/s/F98s5WEo9xMPod7/download"),
      Samoylov_ERA_MkL3_CCSM4_long_term = Resource("Samoylov_ERA_MkL3_CCSM4_long_term", "json", "https://nextcloud.awi.de/s/RSqBtp5sPwkCf45/download"),
)

"""
      Forcing{T,N}

Abstract type representing a generic external boundary condition (i.e. "forcing").
"""
abstract type Forcing{T,N} end
(forcing::Forcing)(x::Number) = error("$(typeof(forcing)) not implemented")
(forcing::Forcing)(t::DateTime) = forcing(ustrip(u"s", float(Dates.datetime2epochms(t))u"ms"))

"""
      Forcings{F<:NamedTuple}

Convenience container for forcings that prevents forcing fields from being included in automatic
type flattening/reconstruction.
"""
@flattenable struct Forcings{F<:NamedTuple}
      forcings::F | false
      Forcings(forcings::F) where {F<:NamedTuple} = new{F}(forcings)
      function Forcings(;kwargs...)
            forcings = (;kwargs...)
            new{typeof(forcings)}(forcings)
      end
end
Base.getproperty(f::Forcings, name::Symbol) = Base.getproperty(getfield(f, :forcings), name)

"""
      TimeSeriesForcing{T,A,I}

Forcing provided by a discrete time series of data.
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

"""
Get interpolated forcing value at t seconds from t0.
"""
(forcing::TimeSeriesForcing)(t::Number) = forcing.interp(t) # extract interpolation and evaluate

Base.getindex(f::TimeSeriesForcing, i) = forcing.tarray[i]
function Base.getindex(f::TimeSeriesForcing{T,A,I}, range::StepRange{DateTime,TStep}) where {T,A,I,TStep}
      order(::Interpolations.GriddedInterpolation{T1,N,T2,Gridded{Torder}}) where {T1,N,T2,Torder} = Torder()
      subseries = f.tarray[range]
      TimeSeriesForcing(values(subseries), timestamp(subseries), colnames(f.tarray)[1]; interp=order(f.interp))
end
