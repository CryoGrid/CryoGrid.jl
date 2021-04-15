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
struct SFCC{F,∇F,S} <: FreezeCurve
    f::F # freeze curve function f: (T,...) -> θ
    ∇f::∇F # derivative of freeze curve function
    solver::S # solver for H -> T or T -> H
    SFCC(f::F,∇f::∇F,s::S) where {F<:SFCCFunction,∇F<:Function,S<:SFCCSolver} = new{F,∇F,S}(f,∇f,s)
end

export SFCC

"""
    SFCC(f::SFCCFunction, s::SFCCSolver=SFCCNewtonSolver())

Convenience constructor for SFCC that automatically generates an analytical derivative of the given
freeze curve function `f` using ModelingToolkit/Symbolics.jl. To avoid symbolic tracing issues, the
function should 1) be pure (no side effects or non-mathematical behavior) and 2) avoid indeterminate
control flow such as if-else or while blocks (technically should work but sometimes doesn't...).
Conditional logic can be incorporated via `IfElse.ifelse`. See the documentation for `Symbolics.jl`
for more information and technical details.
"""
function SFCC(f::SFCCFunction, s::SFCCSolver=SFCCNewtonSolver(); dvar=:T, choosefn=first)
    ∇f = generate_derivative(f, dvar; choosefn=choosefn)
    # we wrap ∇f with Base.splat here to avoid a weird issue with in-place splatting causing allocations
    # when applied to runtime generated functions.
    SFCC(f, Base.splat(∇f), s)
end

# Join the declared state variables of the SFCC function and the solver
variables(sfcc::SFCC) = tuplejoin(variables(sfcc.f), variables(sfcc.solver))
"""
Updates T, θ, and C according to the specified SFCC function and solver.
By default, this is implemented as a simple passthrough to the solver.
"""
(sfcc::SFCC)(soil::Soil, heat::Heat, state) = sfcc.solver(soil, heat, state, sfcc.f, sfcc.∇f)

"""
    params(f::SFCCFunction, soil::Soil, heat::Heat, state)

Retrieves a tuple of values corresponding to each parameter declared by SFCCFunction `f` given the
Soil layer, Heat process, and model state. The order of parameters *must match* the argument order
of the freeze curve function `f`.
"""
params(f::SFCCFunction, soil::Soil, heat::Heat, state) = ()
# Fallback implementation of variables for SFCCFunction
variables(f::SFCCFunction) = ()
variables(s::SFCCSolver) = ()

export params

"""
    VanGenuchten <: SFCCFunction

van Genuchten MT, 1980. A closed-form equation for predicting the hydraulic conductivity of unsaturated soils.
    Soil Science Society of America Journal, 44(5): 892–898. DOI: 10.2136/sssaj 1980.03615995004400050002x.

Dall'Amico M, 2010. Coupled water and heat transfer in permafrost modeling. Ph.D. Thesis, University of Trento, pp. 43.
"""
@with_kw struct VanGenuchten <: SFCCFunction
    θres::Float64 = 0.0 # residual water content
    g::Float64 = 9.80665 # acceleration due to gravity
end
variables(::VanGenuchten) = (Parameter(:α), Parameter(:n), Parameter(:Tₘ))
params(f::VanGenuchten, soil::Soil, heat::Heat, state) = (
    state.params.α |> getscalar, 
    state.params.n |> getscalar,
    state.params.Tₘ |> getscalar,
    state.θw,
    state.θp, # θ saturated = porosity
    heat.params.L, # specific latent heat of fusion, L
)
function (f::VanGenuchten)(T,α,n,Tₘ,θtot,θsat,L)
    let θres = f.θres,
        g = f.g,
        m = 1-1/n,
        ψ₀ = 0.0, #(-1/α)*(((θtot-θres)/(θsat-θres))^(-1/m)-1)^(1/n),
        Tstar = Tₘ + g*Tₘ/L*ψ₀,
        ψ(T) = ψ₀ + L/(g*Tstar)*(T-Tstar)*heaviside(Tstar-T); # pressure head at T
        θl = θres + (θsat - θres)*(1 + (-α*ψ(T))^n)^(-m) # van Genuchten
    end
end

export VanGenuchten

"""
    McKenzie <: SFCCFunction

McKenzie JM, Voss CI, Siegel DI, 2007. Groundwater flow with energy transport and water-ice phase change:
    numerical simulations, benchmarks, and application to freezing in peat bogs. Advances in Water Resources,
    30(4): 966–983. DOI: 10.1016/j.advwatres.2006.08.008.
"""
@with_kw struct McKenzie <: SFCCFunction
    θres::Float64 = 0.0 # residual water content
end
variables(::McKenzie) = (Parameter(:δ),)
params(f::McKenzie, soil::Soil, heat::Heat, state) = (
    state.params.γ |> getscalar, 
    state.θw,
    state.θp,
)
function (f::McKenzie)(T,γ,θtot,θsat)
    let θres = f.θres,
        Tref = 273.15; # reference T in K
        # TODO: perhaps T<=Tref should be T<=Tₘ as in VG curve?
        IfElse.ifelse(T<=Tref, θres + (θsat-θres)*exp(-((T - Tref)/γ)^2), θtot)
    end
end

export McKenzie

"""
    Westermann <: SFCCFunction

Westermann, S., Boike, J., Langer, M., Schuler, T. V., and Etzelmüller, B.: Modeling the impact of
    wintertime rain events on the thermal regime of permafrost, The Cryosphere, 5, 945–959,
    https://doi.org/10.5194/tc-5-945-2011, 2011. 
"""
@with_kw struct Westermann <: SFCCFunction
    θres::Float64 = 0.0 # residual water content
end
variables(::Westermann) = (Parameter(:δ),)
params(f::Westermann, soil::Soil, heat::Heat, state) = (
    state.params.δ |> getscalar, 
    state.θw,
)
function (f::Westermann)(T,δ,θtot)
    let θres = f.θres,
        Tref = 273.15; # reference T in K
        # TODO: perhaps T<=Tref should be T<=Tₘ as in VG curve?
        IfElse.ifelse(T<=Tref, θres - (θtot-θres)*(δ/(T-Tref-δ)), θtot)
    end
end

export Westermann

"""
Specialized implementation of Newton's method with backtracking line search for resolving
the energy conservation law, H = TC + Lθ. Attempts to find the root of the corresponding
temperature residual: ϵ = T - (H - Lθ(T)) / C(θ(T)) and uses backtracking to avoid
jumping over the solution. This prevents convergence issues that arise due to discontinuities
and non-monotonic behavior in most common soil freeze curves.
"""
@with_kw struct SFCCNewtonSolver <: SFCCSolver
    maxiter::Int = 10 # maximum number of iterations
    tol::Float64 = 0.01 # absolute tolerance for convergence
    α₀::Float64 = 1.0 # initial step size multiplier
    τ::Float64 = 0.75 # step size decay for backtracking
    onfail::Symbol = Symbol("warn") # error, warn, or ignore
end
convergencefailure(sym::Symbol, i, maxiter, res) = convergencefailure(Val{sym}(), i, maxiter, res)
convergencefailure(::Val{:error}, i, maxiter, res) = error("grid cell $i failed to converge after $maxiter iterations; residual: $(res); You may want to increase 'maxiter' or decrease your integrator step size.")
convergencefailure(::Val{:warn}, i, maxiter, res) = @warn "grid cell $i failed to converge after $maxiter iterations; residual: $(res); You may want to increase 'maxiter' or decrease your integrator step size."
convergencefailure(::Val{:ignore}, i, maxiter, res) = nothing
# Helper function for handling arguments to freeze curve function, f;
# select calls getindex(i) for all array-typed arguments leaves non-array arguments as-is.
# we use a generated function to expand the arguments into an explicitly defined tuple to preserve type-stability (i.e. it's an optmization);
# function f is then applied to each element
@generated selectat(i::Int, f, args::T) where {T<:Tuple} = :(tuple($([typ <: Vector ?  :(f(args[$k][i])) : :(f(args[$k])) for (k,typ) in enumerate(Tuple(T.parameters))]...)))
# Newton solver implementation
function (s::SFCCNewtonSolver)(soil::Soil, heat::Heat{u"J"}, state, f, ∇f)
    # Helper function for updating θl, C, and the residual.
    function residual(T, Tref, H, θw, θm, θo, L, hcparams, f, f_args)
        args = tuplejoin((T,),f_args)
        θl = f(args...)
        C = heatcapacity(hcparams, θw, θl, θm, θo)
        Tres = (T-Tref) - (H - θl*L) / C
        return Tres, θl, C
    end
    # get f arguments; note that this does create some redundancy in the arguments
    # eventually passed to the `residual` function; this is less than ideal but
    # probably shouldn't incur too much of a performance hit, just a few extra stack pointers!
    f_args = params(f, soil, heat, state)
    # iterate over each cell and solve the conservation law: H = TC + Lθ
    @inbounds @fastmath for i in 1:length(state.T)
        itercount = 0
        let Tref = 273.15, # reference temperature: K -> °C
            T₀ = state.T[i] |> adstrip,
            T = T₀, # temperature (K)
            H = state.H[i] |> adstrip, # enthalpy
            C = state.C[i] |> adstrip, # heat capacity
            θl = state.θl[i] |> adstrip, # liquid water content
            θtot = state.θw[i] |> adstrip, # total water content
            θm = state.θm[i] |> adstrip, # mineral content
            θo = state.θo[i] |> adstrip, # organic content
            θp = state.θp[i] |> adstrip, # porosity and/or θsat
            L = heat.params.L, # specific latent heat of fusion
            cw = soil.hcparams.cw, # heat capacity of liquid water
            α₀ = s.α₀,
            τ = s.τ,
            f_argsᵢ = selectat(i, adstrip, f_args);
            # compute initial residual
            Tres, θl, C = residual(T, Tref, H, θtot, θm, θo, L, soil.hcparams, f, f_argsᵢ)
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
                T̂res, θl, C = residual(T̂, Tref, H, θtot, θm, θo, L, soil.hcparams, f, f_argsᵢ)
                # simple backtracking line search to avoid jumping over the solution
                while sign(T̂res) != sign(Tres)
                    α = α*τ # decrease step size by τ
                    T̂ = T - α*Tres # new guess for T
                    T̂res, θl, C = residual(T̂, Tref, H, θtot, θm, θo, L, soil.hcparams, f, f_argsᵢ)
                end
                T = T̂ # update T
                Tres = T̂res # update residual
                itercount += 1
            end
            # Here we apply the optimized result to the state variables;
            # Since we perform the Newton iteration on untracked variables,
            # we need to recompute θl, C, and T here with the tracked variables.
            # Note that this results in one additional freeze curve function evaluation.
            let f_argsᵢ = selectat(i,identity,f_args);
                # recompute liquid water content with (possibly) tracked variables
                args = tuplejoin((T,),f_argsᵢ)
                state.θl[i] = f(args...)
            end
            let θl = state.θl[i],
                H = state.H[i];
                state.C[i] = heatcapacity(soil.hcparams,θtot,θl,θm,θo)
                state.T[i] = (H - L*θl) / state.C[i] + Tref
            end
        end
    end
    nothing
end

export SFCCNewtonSolver
