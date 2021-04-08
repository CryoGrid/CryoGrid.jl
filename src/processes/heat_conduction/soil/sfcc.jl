"""
    SFCC{F,∇F,S} <: FreezeCurve

Generic representation of the soil freeze characteristic curve. The shape and parameters
of the curve are determined by the implementation of SFCCFunction `f`. Also requires
an implementation of SFCCSolver which provides the solution to the non-linear mapping H <--> T.
"""
struct SFCC{F,∇F,S} <: FreezeCurve
    f::F # freeze curve function f: (T,...) -> θ
    ∇f::∇F # derivative of freeze curve function
    solver::S # solver for H -> T or T -> H
    SFCC(f::F,∇f::∇F,s::S) where {F<:SFCCFunction,∇F<:Function,S<:SFCCSolver} = new{F,∇F,S}(f,∇f,solver)
end

"""
Convenience constructor for SFCC that automatically generates an analytical derivative of the given
freeze curve function `f` using ModelingToolkit/Symbolics.jl. To avoid symbolic tracing issues, the
function should 1) be pure (no side effects or non-mathematical behavior) and 2) avoid indeterminate
control flow such as if-else or while blocks (technically should work but sometimes doesn't...).
Conditional logic can be incorporated via `IfElse.ifelse`. See the documentation for `Symbolics.jl`
for more information and technical details.
"""
function SFCC(f::SFCCFunction, s::SFCCSolver)
    # 1. Parse function parameter names using ExprTools
    fms = methods(f)
    @assert length(fms) == 1 "Freeze curve function can only have one method dispatch to avoid ambiguity. Use Base.deletemethod or restart Julia."
    symbol(arg::Symbol) = arg
    symbol(expr::Expr) = expr.args[1]
    argnames = map(symbol, ExprTools.signature(first(fms))[:args])
    @assert length(argnames) >= 1 "Freeze curve function must have at least one argument, T"
    @assert first(argnames) == :T, "Expected $(:T) as the first argument to freeze curve function but found $(first(argnames))"
    # 2. Convert to MTK symbols
    argsyms = map(s -> Num(Sym{Real}(s)), argnames)
    # 3. Generate analytical derivative of f
    @variables T
    ∂T = Differential(T)
    ∇f_expr = build_function(∂T(f(T)) |> expand_derivatives,T)
    ∇f = eval(dFdT_expr)
    SFCC(f, ∇f, s)
end

# Join the declared state variables of the SFCC function and the solver
variables(sfcc::SFCC) = tuplejoin(variables(sfcc.f), variables(sfcc.solver))

"""
Abstract representation of a soil freeze characteristic curve (SFCC) function.
Subtypes should be callable structs that implement the freeze curve and contain
any necessary additional constants or configuration options. User-specified parameters
can either be supplied in the struct or declared as model parameters via the `variables`
method.
"""
abstract type SFCCFunction end
"""
    params(f::SFCCFunction, soil::Soil, heat::Heat, state)

Retrieves a tuple of values corresponding to each parameter declared by SFCCFunction `f` given the
Soil layer, Heat process, and model state. The order of parameters *must match* the argument order
of the freeze curve function `f`.
"""
params(f::SFCCFunction, soil::Soil, heat::Heat, state) = ()
# Fallback implementation of variables for SFCCFunction
variables(f::SFCCFunction) = ()
"""
    VanGenuchtenFunction <: SFCCFunction
"""
@with_kw struct VanGenuchtenFunction <: SFCCFunction
    Tₘ::Float"K" = 273.15 # freezing point of water
    θres::Float64 = 0.0 # residual water content
    g::Float64 = 9.80665 # acceleration due to gravity
end
variables(::VanGenuchtenFunction) = (Parameter(:α), Parameter(:n))
params(f::VanGenuchtenFunction, soil::Soil, heat::Heat, state) = (
    state.T,
    state.params.α, 
    state.params.n,
    state.θw,
    state.θp, # θ saturated = porosity
    heat.params.ρ*heat.params.Lsl, # specific latent heat of fusion, L
)
function (f::VanGenuchtenFunction)(T,α,n,θtot,θsat,L)
    let Tₘ = f.Tₘ,
        θres = f.θres,
        g = f.g,
        m = 1-1/n,
        ψ₀ = (-1/α)*(((θtot-θres)/(θsat-θres))^(-1/m)-1)^(1/n),
        Tstar = Tₘ + g*Tₘ/L*ψ₀,
        ψ(T) = ψ₀ + L/(g*Tstar)*(T-Tstar)*heaviside(Tstar-T); # pressure head at T
        θw = θres + (θsat - θres)*(1 + (-α*ψ(T))^n)^(-m) # van Genuchten
    end
end

"""
    McKenzieFunction <: SFCCFunction
"""
@with_kw struct McKenzieFunction <: SFCCFunction
    θres::Float64 = 0.0 # residual water content
end
variables(::McKenzieFunction) = (Parameter(:δ),)
params(f::McKenzieFunction, soil::Soil, heat::Heat, state) = (state.params.δ, state.θp)
function (f::McKenzieFunction)(T,δ,θsat)
    IfElse.ifelse(T<=0.0, θres + (θsat-θres)*exp(-(T/δ)^2),θtot)
end

"""
Abstract type for SFCC H <--> T solvers.
"""
abstract type SFCCSolver end
"""
Specialized implementation of Newton's method with backtracking line search for realizing
the energy conservation law, H = TC + Lθ. Attempts to find the root of the corresponding
temperature residual: ϵ = T - (H - Lθ(T)) / C(θ(T)) and uses backtracking to avoid
jumping over the solution. This prevents convergence issues that arise due to discontinuities
and non-monotonic behavior in most common soil freeze curves.
"""
@with_kw struct NewtonsMethod <: SFCCSolver
    maxiter::Int = 10 # maximum number of iterations
    tol::float64 = 0.01 # absolute tolerance for convergence
    α₀::Float64 = 1.0 # initial step size multiplier
    τ::Float64 = 0.7 # step size decay for backtracking
end
