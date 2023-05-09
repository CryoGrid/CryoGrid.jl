"""
Analytical solution to linear form of the heat equation on a semi-infinite rod with periodic fluctuations at the upper  boundary.
Solution taken from Riseborough et al. 2008.
"""
function heat_conduction_linear_periodic_ub(T₀, A, P, α)
    T(z,t) = T₀ + A*exp(-z*sqrt(π/(α*P)))*sin(2π*t/P - z*sqrt(π/(α*P)))
    return T
end

