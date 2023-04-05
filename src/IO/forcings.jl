"""
      Forcing{unit,T}

Abstract type representing a generic external boundary condition (i.e. "forcing").
"""
abstract type Forcing{unit,T} end
@inline @propagate_inbounds (forcing::Forcing)(x::Number) = error("$(typeof(forcing)) not implemented")
@inline @propagate_inbounds (forcing::Forcing)(t::DateTime) = forcing(ustrip(u"s", float(Dates.datetime2epochms(t))u"ms"))

struct ConstantForcing{unit,T} <: Forcing{unit,T}
      value::T
      ConstantForcing(qty::Unitful.AbstractQuantity) = new{unit(qty),typeof(qty)}(qty)
      ConstantForcing(qty::Number) = new{Unitful.NoUnits,typeof(qty)}(qty)
end
(f::ConstantForcing)(t::Number) = f.value

"""
      TimeSeriesForcing{unit,T,A,I}

Forcing provided by a discrete time series of data.
"""
struct TimeSeriesForcing{unit,T,A,I} <: Forcing{unit,T}
      tarray::TimeArray{T,1,DateTime,A}
      interpolant::I
      TimeSeriesForcing(tarray::TimeArray{T,N,D,A}, interpolant) where {T,N,D,A} = new{unit(eltype(values(tarray))),T,A,typeof(interpolant)}(tarray, interpolant)
      function TimeSeriesForcing(values::A, timestamps::AbstractArray{DateTime,1}, name::Symbol; interpolation_mode=Linear()) where {T,A<:AbstractArray{T,1}}
            tarray = TimeArray(collect(timestamps), collect(values), [name])
            TimeSeriesForcing(tarray; interpolation_mode)
      end
      function TimeSeriesForcing(tarray::TimeArray{T,N,D,A}; interpolation_mode=Linear()) where {T,N,D,A}
            # converts timestamps to seconds since epoch
            ts = ustrip.(u"s", float.(Dates.datetime2epochms.(timestamp(tarray)))u"ms")
            u = unit(eltype(values(tarray)))
            interpolant = interpolate((ts,), ustrip.(values(tarray)), Gridded(interpolation_mode))
            new{u,T,A,typeof(interpolant)}(tarray, interpolant)
      end
end
CryoGrid.parameterize(f::TimeSeriesForcing; props...) = f
Flatten.flattenable(::Type{<:TimeSeriesForcing}, ::Type) = false

Base.show(io::IO, forcing::TimeSeriesForcing{u}) where u = print(io, "TimeSeriesForcing $(first(colnames(forcing.tarray))) [$u] of length $(length(forcing.tarray)) with time span $(extrema(timestamp(forcing.tarray)))")

"""
Get interpolated forcing value at t seconds from t0.
"""
@inline @propagate_inbounds (forcing::TimeSeriesForcing)(t::Number) = forcing.interpolant(t) # extract interpolation and evaluate

@inline @propagate_inbounds Base.getindex(forcing::TimeSeriesForcing, i) = forcing.tarray[i]
function Base.getindex(f::TimeSeriesForcing{u,T,A,I}, range::StepRange{DateTime,TStep}) where {u,T,A,I,TStep}
      order(::Interpolations.GriddedInterpolation{T1,N,T2,Gridded{Torder}}) where {T1,N,T2,Torder} = Torder()
      subseries = f.tarray[range]
      TimeSeriesForcing(values(subseries), timestamp(subseries), colnames(f.tarray)[1]; interpolation_mode=order(f.interpolant))
end
