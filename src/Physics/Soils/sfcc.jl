"""
Abstract representation of a soil freeze characteristic curve (SFCC) function.
Subtypes should be callable structs that implement the freeze curve and contain
any necessary additional constants or configuration options. User-specified parameters
can either be supplied in the struct or declared as model parameters via the `variables`
method.
"""
abstract type SFCCFunction <: Function end
"""
Abstract type for SFCC H <--> T solvers.
"""
abstract type SFCCSolver end
"""
    SFCC{F,∇F,S} <: FreezeCurve

Generic representation of the soil freeze characteristic curve. The shape and parameters
of the curve are determined by the implementation of SFCCFunction `f`. Also requires
an implementation of SFCCSolver which provides the solution to the non-linear mapping H <--> T.
"""
@flattenable struct SFCC{F,∇F,S} <: FreezeCurve
    f::F | true # freeze curve function f: (T,...) -> θ
    ∇f::∇F | false # derivative of freeze curve function
    solver::S | true # solver for H -> T or T -> H
    SFCC(f::F,∇f::∇F,s::S) where {F<:SFCCFunction,∇F<:Function,S<:SFCCSolver} = new{F,∇F,S}(f,∇f,s)
end

"""
    SFCC(f::SFCCFunction, s::SFCCSolver=SFCCNewtonSolver())

Convenience constructor for SFCC that automatically generates an analytical derivative of the given
freeze curve function `f` using ModelingToolkit/Symbolics.jl. To avoid symbolic tracing issues, the
function should 1) be pure (no side effects or non-mathematical behavior) and 2) avoid indeterminate
control flow such as if-else or while blocks (technically should work but sometimes doesn't...).
Conditional logic can be incorporated via `IfElse.ifelse`. See the documentation for `Symbolics.jl`
for more information and technical details.
"""
function SFCC(f::SFCCFunction, s::SFCCSolver=SFCCNewtonSolver(); dvar=:T, choosefn=first, context_module=Numerics)
    ∇f = ∇(pstrip(f), dvar; choosefn=choosefn, context_module=context_module)
    # we wrap ∇f with Base.splat here to avoid a weird issue with in-place splatting causing allocations
    # when applied to runtime generated functions.
    SFCC(f, ∇f, s)
end

# Join the declared state variables of the SFCC function and the solver
variables(soil::Soil, heat::Heat, sfcc::SFCC) = tuplejoin(variables(soil, heat, sfcc.f), variables(soil, heat, sfcc.solver))
# Default SFCC initialization
function initialcondition!(soil::Soil, heat::Heat, sfcc::SFCC, state)
    L = heat.L
    state.θl .= sfcc.f.(state.T, sfccparams(sfcc.f, soil, heat, state)...)
    heatcapacity!(soil, heat, state)
    @. state.H = enthalpy(state.T, state.C, L, state.θl)
end

"""
    sfccparams(f::SFCCFunction, soil::Soil, heat::Heat, state)

Retrieves a tuple of values corresponding to each parameter declared by SFCCFunction `f` given the
Soil layer, Heat process, and model state. The order of parameters *must match* the argument order
of the freeze curve function `f`.
"""
sfccparams(::SFCCFunction, ::Soil, ::Heat, state) = ()
# Fallback implementation of variables for SFCCFunction
variables(::Soil, ::Heat, f::SFCCFunction) = ()
variables(::Soil, ::Heat, s::SFCCSolver) = ()
"""
    DallAmico <: SFCCFunction

Dall'Amico M, 2010. Coupled water and heat transfer in permafrost modeling. Ph.D. Thesis, University of Trento, pp. 43.
"""
Base.@kwdef struct DallAmico{T,Θ,A,N,G} <: SFCCFunction
    Tₘ::T = Param(0.0, units=u"°C")
    θres::Θ = Param(0.0, bounds=(0,1))
    α::A = Param(4.0, bounds=(eps(),Inf), units=u"1/m")
    n::N = Param(2.0, bounds=(1,Inf))
    g::G = 9.80665u"m/s^2" # acceleration due to gravity
    swrc::VanGenuchten = VanGenuchten()
end
sfccparams(f::DallAmico, soil::Soil, heat::Heat, state) = (
    f.Tₘ,
    f.θres,
    porosity(soil, state), # θ saturated = porosity
    totalwater(soil, state), # total water content
    heat.prop.Lf, # specific latent heat of fusion, L
    f.α,
    f.n,
)
# pressure head at T
@inline ψ(T,Tstar,ψ₀,Lf,g) = ψ₀ + Lf/(g*Tstar)*(T-Tstar)*heaviside(Tstar-T)
@inline function (f::DallAmico)(T,Tₘ,θres,θsat,θtot,Lf,α,n)
    let θsat = max(θtot, θsat),
        g = f.g,
        m = 1-1/n,
        Tₘ = normalize_temperature(Tₘ),
        ψ₀ = IfElse.ifelse(θtot < θsat, -1/α*(((θtot-θres)/(θsat-θres))^(-1/m)-1.0)^(1/n), zero(1/α)),
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
    θres::Θ = Param(0.0, bounds=(0,1))
    γ::Γ = Param(0.1, bounds=(eps(),Inf), units=u"K")
end
sfccparams(f::McKenzie, soil::Soil, heat::Heat, state) = (
    f.Tₘ,
    f.θres,
    porosity(soil, state), # θ saturated = porosity
    totalwater(soil, state), # total water content
    f.γ,
)
function (f::McKenzie)(T,Tₘ,θres,θsat,θtot,γ)
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
    θres::Θ = Param(0.0, bounds=(0,1))
    δ::Δ = Param(0.1, bounds=(eps(),Inf), units=u"K")
end
sfccparams(f::Westermann, soil::Soil, heat::Heat, state) = (
    f.Tₘ,
    f.θres,
    porosity(soil, state), # θ saturated = porosity
    totalwater(soil, state), # total water content
    f.δ,
)
function (f::Westermann)(T,Tₘ,θres,θsat,θtot,δ)
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
"""
    SFCC(f::SFCCTable, s::SFCCSolver=SFCCNewtonSolver())

Constructs a SFCC from the precomputed `SFCCTable`. The derivative is generated using the
`gradient` function provided by `Interpolations`.
"""
function SFCC(f::SFCCTable, s::SFCCSolver=SFCCNewtonSolver())
    # we wrap ∇f with Base.splat here to avoid a weird issue with in-place splatting causing allocations
    # when applied to runtime generated functions.
    SFCC(f, Base.splat(first ∘ ∇(f.f_tab)), s)
end

include("sfcc_solvers.jl")

# Generate analytical derivatives during precompilation
const ∂DallAmico∂T = ∇(stripunits(stripparams(DallAmico())), :T)
const ∂McKenzie∂T = ∇(stripunits(stripparams(McKenzie())), :T)
const ∂Westermann∂T = ∇(stripunits(stripparams(Westermann())), :T)
SFCC(f::DallAmico, solver::SFCCSolver=SFCCNewtonSolver()) = SFCC(f, ∂DallAmico∂T, solver)
SFCC(f::McKenzie, solver::SFCCSolver=SFCCNewtonSolver()) = SFCC(f, ∂McKenzie∂T, solver)
SFCC(f::Westermann, solver::SFCCSolver=SFCCNewtonSolver()) = SFCC(f, ∂Westermann∂T, solver)
