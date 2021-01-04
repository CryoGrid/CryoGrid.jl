using DiffEqOperators

const Δ = CenteredDifference(ord_deriv, ord_approx, h*ones(nknots+1), nknots)
const bc = Dirichlet0BC(Float64)

t0 = 0.0
t1 = 1.0
u0 = 200.0*Vector(1.0.-knots./2).-1.0

function dTdt(dT,T,p,t)
    dT .= Δ*bc*T
end
