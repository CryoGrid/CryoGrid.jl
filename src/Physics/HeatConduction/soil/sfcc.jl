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
            state.θl[i] = f(f_argsᵢ...)
            state.C[i] = heatcapacity(soil, state.θw[i], state.θl[i], state.θm[i], state.θo[i])
            state.H[i] = enthalpy(state.T[i], state.C[i], L, state.θl[i])
            state.Ceff[i] = L*∇f(f_argsᵢ) + state.C[i]
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
    porosity(soil), # θ saturated = porosity
    state.θw, # total water content
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
    porosity(soil), # θ saturated = porosity
    state.θw, # total water content
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
    porosity(soil), # θ saturated = porosity
    state.θw, # total water content
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
Numerics.Tabulated(f::SFCCFunction, args...) = SFCCTable(f, Numerics.tabulate(f, args...))
"""
    SFCC(f::SFCCTable, s::SFCCSolver=SFCCNewtonSolver())

Constructs a SFCC from the precomputed `SFCCTable`. The derivative is generated using the
`gradient` function provided by `Interpolations`.
"""
function SFCC(f::SFCCTable, s::SFCCSolver=SFCCNewtonSolver())
    # we wrap ∇f with Base.splat here to avoid a weird issue with in-place splatting causing allocations
    # when applied to runtime generated functions.
    SFCC(f, Base.splat(∇(f.f_tab)), s)
end

"""
Specialized implementation of Newton's method with backtracking line search for resolving
the energy conservation law, H = TC + Lθ. Attempts to find the root of the corresponding
temperature residual: ϵ = T - (H - Lθ(T)) / C(θ(T)) and uses backtracking to avoid
jumping over the solution. This prevents convergence issues that arise due to discontinuities
and non-monotonic behavior in most common soil freeze curves.
"""
@with_kw struct SFCCNewtonSolver <: SFCCSolver
    maxiter::Int = 50 # maximum number of iterations
    tol::Float64 = 0.01 # absolute tolerance for convergence
    α₀::Float64 = 1.0 # initial step size multiplier
    τ::Float64 = 0.7 # step size decay for backtracking
    onfail::Symbol = Symbol("warn") # error, warn, or ignore
end
convergencefailure(sym::Symbol, i, maxiter, res) = convergencefailure(Val{sym}(), i, maxiter, res)
convergencefailure(::Val{:error}, i, maxiter, res) = error("grid cell $i failed to converge after $maxiter iterations; residual: $(res); You may want to increase 'maxiter' or decrease your integrator step size.")
convergencefailure(::Val{:warn}, i, maxiter, res) = @warn "grid cell $i failed to converge after $maxiter iterations; residual: $(res); You may want to increase 'maxiter' or decrease your integrator step size."
convergencefailure(::Val{:ignore}, i, maxiter, res) = nothing
# Helper function for updating θl, C, and the residual.
function residual(T, H, θw, θm, θo, L, soil, f, f_args)
    args = tuplejoin((T,),f_args)
    θl = f(args...)
    C = heatcapacity(soil, θw, θl, θm, θo)
    Tres = T - (H - θl*L) / C
    return Tres, θl, C
end
# Newton solver implementation
function (s::SFCCNewtonSolver)(soil::Soil, heat::Heat{<:SFCC,Enthalpy}, state, f, ∇f)
    # get f arguments; note that this does create some redundancy in the arguments
    # eventually passed to the `residual` function; this is less than ideal but
    # probably shouldn't incur too much of a performance hit, just a few extra stack pointers!
    f_args = sfccparams(f, soil, heat, state)
    # iterate over each cell and solve the conservation law: H = TC + Lθ
    @threaded for i in 1:length(state.T)
        @inbounds @fastmath let T₀ = state.T[i] |> Utils.adstrip,
            T = T₀, # temperature
            H = state.H[i] |> Utils.adstrip, # enthalpy
            C = state.C[i] |> Utils.adstrip, # heat capacity
            θl = state.θl[i] |> Utils.adstrip, # liquid water content
            θw = state.θw[i] |> Utils.adstrip, # total water content
            θm = state.θm[i] |> Utils.adstrip, # mineral content
            θo = state.θo[i] |> Utils.adstrip, # organic content
            L = heat.L, # specific latent heat of fusion
            cw = soil.hc.cw, # heat capacity of liquid water
            α₀ = s.α₀,
            τ = s.τ,
            f_argsᵢ = Utils.selectat(i, Utils.adstrip, f_args),
            itercount = 0;
            # compute initial guess T by setting θl according to free water scheme
            T = let Lθ = L*θw;
                if H < 0
                    H / heatcapacity(soil,θw,0.0,θm,θo)
                elseif H >= 0 && H < Lθ
                    (1.0 - H/Lθ)*0.1
                else
                    (H - Lθ) / heatcapacity(soil,θw,θw,θm,θo)
                end
            end
            # compute initial residual
            Tres, θl, C = residual(T, H, θw, θm, θo, L, soil, f, f_argsᵢ)
            while abs(Tres) > s.tol
                if itercount > s.maxiter
                    convergencefailure(s.onfail, i, s.maxiter, Tres)
                    break
                end
                # derivative of freeze curve
                args = tuplejoin((T,),f_argsᵢ)
                ∂θ∂T = ∇f(args)
                # derivative of residual by quotient rule;
                # note that this assumes heatcapacity to be a simple weighted average!
                # in the future, it might be a good idea to compute an automatic derivative
                # of heatcapacity in addition to the freeze curve function.
                ∂Tres∂T = 1.0 - ∂θ∂T*(-L*C - (H - θl*L)*cw)/C^2
                α = α₀ / ∂Tres∂T
                T̂ = T - α*Tres
                # do first residual check outside of loop;
                # this way, we don't decrease α unless we have to.
                T̂res, θl, C = residual(T̂, H, θw, θm, θo, L, soil, f, f_argsᵢ)
                inneritercount = 0
                # simple backtracking line search to avoid jumping over the solution
                while sign(T̂res) != sign(Tres)
                    if inneritercount > 100
                        @warn "Backtracking failed; this should not happen. Current state: α=$α, T=$T, T̂=$T̂, residual $(T̂res), initial residual: $(Tres)"
                        break
                    end
                    α = α*τ # decrease step size by τ
                    T̂ = T - α*Tres # new guess for T
                    T̂res, θl, C = residual(T̂, H, θw, θm, θo, L, soil, f, f_argsᵢ)
                    inneritercount += 1
                end
                T = T̂ # update T
                Tres = T̂res # update residual
                itercount += 1
            end
            # Here we apply the optimized result to the state variables;
            # Since we perform the Newton iteration on untracked variables,
            # we need to recompute θl, C, and T here with the tracked variables.
            # Note that this results in one additional freeze curve function evaluation.
            let f_argsᵢ = Utils.selectat(i,identity,f_args);
                # recompute liquid water content with (possibly) tracked variables
                args = tuplejoin((T,),f_argsᵢ)
                state.θl[i] = f(args...)
                dθdT = ∇f(args)
                let θl = state.θl[i],
                    H = state.H[i];
                    state.C[i] = heatcapacity(soil,θw,θl,θm,θo)
                    state.Ceff[i] = state.C[i] + dθdT
                    state.T[i] = (H - L*θl) / state.C[i]
                end
            end
        end
    end
    nothing
end

# Generate analytical derivatives during precompilation
const ∂DallAmico∂T = ∇(DallAmico(), :T)
const ∂McKenzie∂T = ∇(McKenzie(), :T)
const ∂Westermann∂T = ∇(Westermann(), :T)
SFCC(f::DallAmico, solver::SFCCSolver=SFCCNewtonSolver()) = SFCC(f, Base.splat(∂DallAmico∂T), solver)
SFCC(f::McKenzie, solver::SFCCSolver=SFCCNewtonSolver()) = SFCC(f, Base.splat(∂McKenzie∂T), solver)
SFCC(f::Westermann, solver::SFCCSolver=SFCCNewtonSolver()) = SFCC(f, Base.splat(∂Westermann∂T), solver)
