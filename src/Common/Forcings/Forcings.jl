module Forcings

using Dates
using Interpolations
using TimeSeries
using Unitful

export Forcing, TimeSeriesForcing

"""
      Forcing{T,N}
"""
abstract type Forcing{T,N} end

include("timeseries.jl")

end
