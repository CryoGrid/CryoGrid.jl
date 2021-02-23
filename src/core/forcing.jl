struct Forcing{T,A}
      tarray::TimeArray{T,1,DateTime,A}
      Forcing(values::A, timestamps::AbstractArray{DateTime,1},
              name::Symbol=:forcing; interpolation=Linear) where
              {T,A<:AbstractArray{T,1}} = begin
              degree = interpolation()
              ts = datetime2unix.(timestamps)
              interp = interpolate((ts,), values, Gridded(degree))
              new{T,typeof(interp)}(TimeArray(timestamps, interp, [name]))
      end
end

Base.getindex(f::Forcing, t::DateTime) = f[datetime2unix(t)]
Base.getindex(f::Forcing, t::Float64) = values(f.tarray)(t)
Base.getindex(f::Forcing, ts::AbstractArray{DateTime,1}) = getindex.(f,datetime2unix.(ts))
Base.getindex(f::Forcing, ts::AbstractArray{Float64,1}) = values(f.tarray).(ts)
