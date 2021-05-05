########
# Simple heat conduction model (no forcings)
# Uses ModelingToolkit
# Currently does not work because the MOL discretization strategy does not
# support variables without time derivatives, which is required for this model.
########

using DifferentialEquations
using DiffEqOperators
using ModelingToolkit
using IfElse

ivars = @parameters t z L θtot θm θo
dvars = @variables T(..) H(..) Hs(..) Hl(..) θ(..) C(..) k(..)
Dz = Differential(z)
Dt = Differential(t)
heatcap(θ,θtot,θm,θo) = heatcapacity(SoilHCParams(), θtot, θ, θm, θo)
thermalcond(θ,θtot,θm,θo) = thermalconductivity(SoilTCParams(), θtot, θ, θm, θo)
# free water freeze curve
dθdT(T,L,θtot) = IfElse.ifelse(T >= 0.0, IfElse.ifelse(T <= L*θtot, 1e8, 0.0), 0.0)
@register heatcap(θ,θtot,θm,θo)
@register thermalcond(θ,θtot,θm,θo)
@register dθdT(T,L,θtot)

function mtk_hc_phase_change_simple(dz = range(0.0,1000.0,length=100))
    eqs = [
        θ(t,z) ~ Hl / L,
        C(t,z) ~ heatcap(θ,θtot,θm,θo),
        T(t,z) ~ Hs / C,
        H(t,z) ~ Hs + Hl,
        k(t,z) ~ thermalcond(θ,θtot,θm,θo),
        Dt(H) ~ Dz(k(t,z)*Dz(T(t,z))),
        Dt(Hs) ~ Dt(H(t,z)) / (L/C*dθdT(T,L,θtot) + 1),
        Dt(Hl) ~ Dt(H) - Dt(Hs),
    ]
    # Space and time domains
    domains = [t ∈ IntervalDomain(0.0,24*3600.0),
               z ∈ IntervalDomain(0.0,1000.0)]
    bcs = [
        Hs(z,0.0) ~ -2e6,
        Hs(0.0,t) ~ 2e6,
        Hs(1000.0,t) ~ 0.05,
        Hl(z,0.0) ~ 0.0,
        Hl(0.0,t) ~ 0.0,
        Hl(1000.0,t) ~ 0.0,
        # T(z,0) ~ sin(z*π),
        # T(0.0,t) ~ 10*sin(2π*t),
        # T(1000.0,t) ~ 10.2
    ]
    pdesys = PDESystem(eqs, bcs, domains, ivars, dvars)
    discretization = MOLFiniteDifference([z=>dz],t)
    return discretize(pdesys, discretization), pdesys
end

function mtk_hc_simple(dz = range(0.0,1000.0,length=100))
    @parameters t z
    @variables T(..) k(..)
    Dz = Differential(z)
    Dt = Differential(t)
    eqs = [
        Dt(T) ~ Dz(k(t,z)*Dz(T(t,z)))
    ]
    # Space and time domains
    domains = [t ∈ IntervalDomain(0.0,24*3600.0),
               z ∈ IntervalDomain(0.0,1000.0)]
    bcs = [
        T(0,z) ~ sin(π*z),
        T(t,0.0) ~ 0.0,
        T(t,1000.0) ~ 10.2,
        k(0,z) ~ 1 + 0.5*sin(z)
    ]
    pdesys = PDESystem(eqs, bcs, domains, [t,z], [T,k])
    discretization = MOLFiniteDifference([z=>dz],t)
    return discretize(pdesys, discretization), pdesys
end
