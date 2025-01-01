"""
    relative_to_specific_humidity(r_h, pr, Ts, Tair)

Derives specific humidity from measured relative humidity, air pressure, and soil/air temperatures.
"""
relative_to_specific_humidity(r_h, pr, Ts, Tair) = 0.622*(r_h/100)*vapor_pressure(Tair, Ts) / pr

# saturation vapor pressure
"""
    vapor_pressure(T, a₁, a₂, a₃)

Saturation vapor pressure as a function of air temperature with empirical
coefficients a₁, a₂, and a₃.
"""
vapor_pressure(T, a₁, a₂, a₃) = a₁*exp(a₂*T/(T+a₃))

# saturation vapor pressure for frozen vs. unfrozen conditions
"""
    vapor_pressure(Tair, Ts)

Saturation vapor pressure from air and soil temperature, accounting for both frozen (<0°C) and unfrozen conditions.
"""
vapor_pressure(Tair, Ts) = Ts < zero(Ts) ? vapor_pressure(Tair, 611.0, 22.46, 272.62) : vapor_pressure(Tair, 611.0, 17.62, 243.12)
