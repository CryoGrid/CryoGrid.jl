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
    ∇f = ∇(f, dvar; choosefn=choosefn, context_module=context_module)
    # we wrap ∇f with Base.splat here to avoid a weird issue with in-place splatting causing allocations
    # when applied to runtime generated functions.
    SFCC(f, Base.splat(∇f), s)
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
Updates state variables according to the specified SFCC function and solver.
For heat conduction with enthalpy, this is implemented as a simple passthrough to the non-linear solver.
For heat conduction with temperature, we can simply evaluate the freeze curve to get C_eff, θl, and H.
"""
(sfcc::SFCC)(soil::Soil, heat::Heat{<:SFCC,Enthalpy}, state) = sfcc.solver(soil, heat, state, sfcc.f, sfcc.∇f)
function (sfcc::SFCC)(soil::Soil, heat::Heat{<:SFCC,Temperature}, state)
    @inbounds @fastmath let L = heat.L,
        f = sfcc.f,
        ∇f = sfcc.∇f,
        f_args = tuplejoin((state.T,),sfccparams(f,soil,heat,state));
        for i in 1:length(state.T)
            f_argsᵢ = Utils.selectat(i, identity, f_args)
            θw = totalwater(soil, heat, state, i)
            θm = mineral(soil, heat, i)
            θo = organic(soil, heat, i)
            state.θl[i] = f(f_argsᵢ...)
            state.C[i] = heatcapacity(soil, θw, state.θl[i], θm, θo)
            state.H[i] = enthalpy(state.T[i], state.C[i], L, state.θl[i])
            state.dHdT[i] = L*∇f(f_argsᵢ) + state.C[i]
        end
    end
    return nothing
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
@with_kw struct DallAmico{T,Θ,A,N} <: SFCCFunction
    Tₘ::T = Param(0.0)
    θres::Θ = Param(0.0, bounds=(0,1))
    α::A = Param(4.0, bounds=(eps(),Inf))
    n::N = Param(2.0, bounds=(1,Inf))
    swrc::VanGenuchten = VanGenuchten()
end
sfccparams(f::DallAmico, soil::Soil, heat::Heat, state) = (
    f.Tₘ,
    f.θres,
    porosity(soil, heat, state), # θ saturated = porosity
    totalwater(soil, heat, state), # total water content
    heat.L, # specific latent heat of fusion, L
    f.α,
    f.n,
)
# pressure head at T
ψ(T,Tstar,ψ₀,L,g) = ψ₀ + L/(g*Tstar)*(T-Tstar)*heaviside(Tstar-T)
function (f::DallAmico)(T,Tₘ,θres,θsat,θtot,L,α,n)
    let θsat = max(θtot, θsat),
        g = 9.80665, # acceleration due to gravity
        m = 1-1/n,
        Tₘ = Tₘ + 273.15,
        ψ₀ = IfElse.ifelse(θtot < θsat, -1/α*(((θtot-θres)/(θsat-θres))^(-1/m)-1)^(1/n), 0),
        T = T + 273.15,
        Tstar = Tₘ + g*Tₘ/L*ψ₀,
        ψ = ψ(T,Tstar,ψ₀,L,g);
        f.swrc(ψ,θres,θsat,α,n)
    end
end

"""
    McKenzie <: SFCCFunction

McKenzie JM, Voss CI, Siegel DI, 2007. Groundwater flow with energy transport and water-ice phase change:
    numerical simulations, benchmarks, and application to freezing in peat bogs. Advances in Water Resources,
    30(4): 966–983. DOI: 10.1016/j.advwatres.2006.08.008.
"""
@with_kw struct McKenzie{T,Θ,Γ} <: SFCCFunction
    Tₘ::T = Param(0.0)
    θres::Θ = Param(0.0, bounds=(0,1))
    γ::Γ = Param(0.1, bounds=(eps(),Inf))
end
sfccparams(f::McKenzie, soil::Soil, heat::Heat, state) = (
    f.Tₘ,
    f.θres,
    porosity(soil, heat, state), # θ saturated = porosity
    totalwater(soil, heat, state), # total water content
    f.γ,
)
function (f::McKenzie)(T,Tₘ,θres,θsat,θtot,γ)
    let θsat = max(θtot, θsat);
        IfElse.ifelse(T<=Tₘ, θres + (θsat-θres)*exp(-(T/γ)^2), θtot)
    end
end

"""
    Westermann <: SFCCFunction

Westermann, S., Boike, J., Langer, M., Schuler, T. V., and Etzelmüller, B.: Modeling the impact of
    wintertime rain events on the thermal regime of permafrost, The Cryosphere, 5, 945–959,
    https://doi.org/10.5194/tc-5-945-2011, 2011. 
"""
@with_kw struct Westermann{T,Θ,Δ} <: SFCCFunction
    Tₘ::T = Param(0.0)
    θres::Θ = Param(0.0, bounds=(0,1))
    δ::Δ = Param(0.1, bounds=(eps(),Inf))
end
sfccparams(f::Westermann, soil::Soil, heat::Heat, state) = (
    f.Tₘ,
    f.θres,
    porosity(soil, heat, state), # θ saturated = porosity
    totalwater(soil, heat, state), # total water content
    f.δ,
)
function (f::Westermann)(T,Tₘ,θres,θsat,θtot,δ)
    let θsat = max(θtot, θsat);
        IfElse.ifelse(T<=Tₘ, θres - (θsat-θres)*(δ/(T-δ)), θtot)
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
const ∂DallAmico∂T = ∇(DallAmico(), :T)
const ∂McKenzie∂T = ∇(McKenzie(), :T)
const ∂Westermann∂T = ∇(Westermann(), :T)
SFCC(f::DallAmico, solver::SFCCSolver=SFCCNewtonSolver()) = SFCC(f, Base.splat(∂DallAmico∂T), solver)
SFCC(f::McKenzie, solver::SFCCSolver=SFCCNewtonSolver()) = SFCC(f, Base.splat(∂McKenzie∂T), solver)
SFCC(f::Westermann, solver::SFCCSolver=SFCCNewtonSolver()) = SFCC(f, Base.splat(∂Westermann∂T), solver)