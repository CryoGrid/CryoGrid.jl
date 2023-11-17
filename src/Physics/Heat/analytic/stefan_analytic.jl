using ComponentArrays
using SpecialFunctions
using UnPack

import SimpleNonlinearSolve
import SciMLBase

const DefaultHeatProperties = HeatBalanceProperties()
const DefaultThermalProperties = ThermalProperties()

Utils.@properties StefanParameters(
    c_l = DefaultThermalProperties.ch_w,
    c_s = DefaultThermalProperties.ch_i,
    k_l = DefaultThermalProperties.kh_w,
    k_s = DefaultThermalProperties.kh_i,
    ρ = DefaultHeatProperties.ρw,
    Lf = DefaultHeatProperties.Lsl,
    θwi = 1.0,
    T_m = 0.0u"°C",
    T_s = 0.0u"°C",
    T_l = 1.0u"°C",
)

# Solution functions according to Neumann's method (Neumann, 1912)
# Reproduced from Hu et al. 1996
"""
    stefan_boundary(x0, t0, λ, α, t)

Evaluates the position of the Stefan moving boundary at time `t` given the Stefan constant `λ`,
diffusivity `α`, and initial time and position `t0`, `x0`.
"""
stefan_boundary(x0, t0, λ, α, t) = x0 + 2*λ*sqrt(α*(t-t0))
"""
    stefan_temperature_liquid(T_l, T_m, λ, α_l, x0, t0, x, t)

Evaluates the generic Stefan analytical solution for temperature in the liquid region.
"""
stefan_temperature_liquid(T_l, T_m, λ, α_l, x0, t0, x, t) = T_l - (T_l - T_m)*erf((x-x0)/2*sqrt(α_l*(t-t0))^-1) / erf(λ)
"""
    stefan_temperature_solid(T_s, T_m, λ, α_s, α_l, x0, t0, x, t)

Evaluates the generic Stefan analytical solution for temperature in the solid region.
"""
stefan_temperature_solid(T_s, T_m, λ, α_s, α_l, x0, t0, x, t) = T_s + (T_m - T_s)*erfc((x-x0)/2*sqrt(α_s*(t-t0))^-1) / erfc(λ*sqrt(α_l/α_s))
"""
    stefan_number(c, ρ, Lf, ΔT)

Calculates the Stefan number for heat capacity `c`, desnity `ρ`, specific latent heat of fusion `Lf`, and temperature delta `ΔT`.
"""
stefan_number(c, ρ, Lf, θwi, ΔT) = c*ΔT/(θwi*Lf*ρ)
"""
    stefan_residual(λ, T_m, T_l, T_s, k_l, c_l, k_s, c_s, ρ, Lf)

Evaluates the residual function for the transcendental equation in the two-phase Stefan problem.
"""
function stefan_residual(λ, T_m, T_l, T_s, k_l, c_l, k_s, c_s, ρ, θwi, Lf)
    let α_l = k_l / c_l,
        α_s = k_s / c_s,
        St_l = stefan_number(c_l, ρ, Lf, θwi, T_l - T_m),
        St_s = stefan_number(c_s, ρ, Lf, θwi, T_m - T_s); # zero for one-phase problem
        # Eq. 13 in Hu et al. 1996
        St_l / (exp(λ^2)*erf(λ)) - St_s*sqrt(α_s) / (sqrt(α_l)*exp(α_l*λ^2/α_s)*erfc(λ*sqrt(α_l/α_s))) - λ*sqrt(π)
    end
end
stefan_one_phase_residual(λ, T_m, T_l, k_l, c_l, ρ, θwi, Lf) = stefan_residual(λ, T_m, T_l, T_m, k_l, c_l, k_l, c_l, ρ, θwi, Lf)

"""
    StefanProblem{Tp<:StefanParameters,Tx,Tt}

Represents the simple two-phase Stefan problem defined on a semi-infinite slab. The one-phase Stefan problem can be
computed by setting the parameters `T_s = T_m`.
"""
Base.@kwdef struct StefanProblem{Tp<:StefanParameters,Tx,Tt}
    p::Tp = StefanParameters()
    x0::Tx = 0.0u"m"
    t0::Tt = 0.0u"s"
end
struct StefanSolution{Tprob,Tsol,Tλ}
    prob::Tprob
    nlsol::Tsol
    λ::Tλ
end
"""
    (sol::StefanSolution)(t)

Moving boundary solution, `x_m(t)`.
"""
(sol::StefanSolution)(t) = stefan_boundary(sol.prob.x0, sol.prob.t0, sol.λ, sol.prob.p.k_l / sol.prob.p.c_l, t)
"""
    (sol::StefanSolution)(x,t)

Temperature solution, `T(x,t)`.
"""
function (sol::StefanSolution)(x, t)
    @unpack x0, t0 = sol.prob
    @unpack T_m, T_l, T_s, k_l, c_l, k_s, c_s = sol.prob.p
    let λ = sol.λ,
        α_l = k_l / c_l,
        α_s = k_s / c_s,
        x_m = stefan_boundary(x0, t0, λ, α_l, t);
        T_l = stefan_temperature_liquid(T_l, T_m, λ, α_l, x0, t0, x, t)
        T_s = stefan_temperature_solid(T_s, T_m, λ, α_s, α_l, x0, t0, x, t)
        IfElse.ifelse(x >= x_m, uconvert(u"°C", T_s), uconvert(u"°C", T_l))
    end
end
function SciMLBase.solve(prob::StefanProblem, alg=SimpleNonlinearSolve.SimpleNewtonRaphson(); p=prob.p, x0=prob.x0, t0=prob.t0)
    prob = StefanProblem(ComponentVector(p), x0, t0)
    f(u,p) = stefan_residual(u, p.T_m, p.T_l, p.T_s, p.k_l, p.c_l, p.k_s, p.c_s, p.ρ, p.θwi, p.Lf)
    nlprob = SimpleNonlinearSolve.NonlinearProblem(f, ustrip(1/(p.Lf*p.θwi)), p)
    nlsol = SciMLBase.solve(nlprob, alg)
    λ = nlsol.u
    return StefanSolution(prob, nlsol, λ)
end
