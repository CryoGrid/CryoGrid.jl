module Forcings

using TimeSeries

export Forcing
export TotalPrecipitation, AirTemperature, Snowfall

abstract type Forcing{T,D<:TimeType} end

function Timestamps(forcing::Forcing{T,D}) where {T,D<:TimeType}
      notimplemented(typeof(forcing))
end

end
