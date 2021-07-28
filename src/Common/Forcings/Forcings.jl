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
(forcing::Forcing)(x::Number) = error("$(typeof(forcing)) not implemented")
(forcing::Forcing)(t::DateTime) = forcing(ustrip(u"s", float(Dates.datetime2epochms(t))u"ms"))

include("timeseries.jl")

end
