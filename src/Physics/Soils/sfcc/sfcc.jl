"""
Abstract representation of a soil freeze characteristic curve (SFCC) function.
Subtypes should be callable structs that implement the freeze curve and contain
any necessary additional constants or configuration options. User-specified parameters
can either be supplied in the struct or declared as model parameters via the `variables`
method.
"""
abstract type SFCCFunction end
"""
Abstract type for SFCC H <--> T solvers.
"""
abstract type SFCCSolver end
"""
    SFCC{F,S} <: FreezeCurve

Generic representation of the soil freeze characteristic curve. The shape and parameters
of the curve are determined by the implementation of SFCCFunction `f`. Also requires
an implementation of SFCCSolver which provides the solution to the non-linear mapping H <--> T.
"""
struct SFCC{F,S} <: FreezeCurve
    f::F # freeze curve function f: (T,...) -> θ
    solver::S # solver for H -> T or T -> H
    SFCC(f::F, s::S=SFCCPreSolver()) where {F<:SFCCFunction,S<:SFCCSolver} = new{F,S}(f,s)
end

# Join the declared state variables of the SFCC function and the solver
CryoGrid.variables(soil::Soil, heat::Heat, sfcc::SFCC) = tuplejoin(CryoGrid.variables(soil, heat, sfcc.f), CryoGrid.variables(soil, heat, sfcc.solver))
# Default SFCC initialization
function CryoGrid.initialcondition!(soil::Soil, heat::Heat, sfcc::SFCC, state)
    HeatConduction.freezethaw!(soil, heat, state)
end

"""
    sfccargs(f::SFCCFunction, soil::Soil, heat::Heat, state)

Retrieves a tuple of values corresponding to each parameter declared by SFCCFunction `f` given the
Soil layer, Heat process, and model state. The order of parameters *must match* the argument order
of the freeze curve function `f`.
"""
sfccargs(::SFCCFunction, ::Soil, ::Heat, state) = ()
# Fallback implementation of variables for SFCCFunction
CryoGrid.variables(::Soil, ::Heat, f::SFCCFunction) = ()
CryoGrid.variables(::Soil, ::Heat, s::SFCCSolver) = ()
"""
    DallAmico <: SFCCFunction

Dall'Amico M, 2010. Coupled water and heat transfer in permafrost modeling. Ph.D. Thesis, University of Trento, pp. 43.
"""
Base.@kwdef struct DallAmico{T,Θ,G,Tvg<:VanGenuchten} <: SFCCFunction
    Tₘ::T = Param(0.0, units=u"°C")
    θres::Θ = Param(0.0, domain=0..1)
    g::G = 9.80665u"m/s^2" # acceleration due to gravity
    swrc::Tvg = VanGenuchten()
end
sfccargs(f::DallAmico, soil::Soil, heat::Heat, state) = (
    porosity(soil, state), # θ saturated = porosity
    waterice(soil, heat, state), # total water content
    heat.prop.Lf, # specific latent heat of fusion, L
    f.θres,
    f.Tₘ,
    f.swrc.α,
    f.swrc.n,
)
# pressure head at T
@inline ψ(T,Tstar,ψ₀,Lf,g) = ψ₀ + Lf/(g*Tstar)*(T-Tstar)*heaviside(Tstar-T)
@inline function (f::DallAmico)(T, θsat, θtot, Lf, θres=f.θres, Tₘ=f.Tₘ, α=f.swrc.α, n=f.swrc.n)
    let θsat = max(θtot, θsat),
        g = f.g,
        Tₘ = normalize_temperature(Tₘ),
        ψ₀ = f.swrc(inv, θtot, θres, θsat, α, n),
        Tstar = Tₘ + g*Tₘ/Lf*ψ₀,
        T = normalize_temperature(T),
        ψ = ψ(T, Tstar, ψ₀, Lf, g);
        f.swrc(ψ, θres, θsat, α, n)
    end
end

"""
    McKenzie <: SFCCFunction

McKenzie JM, Voss CI, Siegel DI, 2007. Groundwater flow with energy transport and water-ice phase change:
    numerical simulations, benchmarks, and application to freezing in peat bogs. Advances in Water Resources,
    30(4): 966–983. DOI: 10.1016/j.advwatres.2006.08.008.
"""
Base.@kwdef struct McKenzie{T,Θ,Γ} <: SFCCFunction
    Tₘ::T = Param(0.0, units=u"°C")
    θres::Θ = Param(0.0, domain=0..1)
    γ::Γ = Param(0.1, domain=0..Inf, units=u"K")
end
sfccargs(f::McKenzie, soil::Soil, heat::Heat, state) = (
    porosity(soil, state), # θ saturated = porosity
    waterice(soil, heat, state), # total water content
    f.θres,
    f.Tₘ,
    f.γ,
)
function (f::McKenzie)(T, θsat, θtot, θres=f.θres, Tₘ=f.Tₘ, γ=f.γ)
    let T = normalize_temperature(T),
        Tₘ = normalize_temperature(Tₘ),
        θsat = max(θtot, θsat);
        IfElse.ifelse(T<=Tₘ, θres + (θsat-θres)*exp(-((T-Tₘ)/γ)^2), θtot)
    end
end

"""
    Westermann <: SFCCFunction

Westermann, S., Boike, J., Langer, M., Schuler, T. V., and Etzelmüller, B.: Modeling the impact of
    wintertime rain events on the thermal regime of permafrost, The Cryosphere, 5, 945–959,
    https://doi.org/10.5194/tc-5-945-2011, 2011. 
"""
Base.@kwdef struct Westermann{T,Θ,Δ} <: SFCCFunction
    Tₘ::T = Param(0.0, units=u"°C")
    θres::Θ = Param(0.0, domain=0..1)
    δ::Δ = Param(0.1, domain=0..Inf, units=u"K")
end
sfccargs(f::Westermann, soil::Soil, heat::Heat, state) = (
    porosity(soil, state), # θ saturated = porosity
    waterice(soil, heat, state), # total water content
    f.θres,
    f.Tₘ,
    f.δ,
)
function (f::Westermann)(T,θsat,θtot,θres=f.θres,Tₘ=f.Tₘ,δ=f.δ)
    let T = normalize_temperature(T),
        Tₘ = normalize_temperature(Tₘ),
        θsat = max(θtot, θsat);
        IfElse.ifelse(T<=Tₘ, θres - (θsat-θres)*(δ/(T-Tₘ-δ)), θtot)
    end
end
struct SFCCTable{F,I} <: SFCCFunction
    f::F
    f_tab::I
end
(f::SFCCTable)(args...) = f.f_tab(args...)
"""
    Tabulated(f::SFCCFunction, args...)

Produces an `SFCCTable` function which is a tabulation of `f`.
"""
Numerics.Tabulated(f::SFCCFunction, args...; kwargs...) = SFCCTable(f, Numerics.tabulate(f, args...; kwargs...))

include("sfcc_solvers.jl")
