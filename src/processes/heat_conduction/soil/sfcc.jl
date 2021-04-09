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
    SFCC(f, ∇f, s)
end

# Join the declared state variables of the SFCC function and the solver
variables(sfcc::SFCC) = tuplejoin(variables(sfcc.f), variables(sfcc.solver))
"""
Updates T, θ, and C according to the specified SFCC function and solver.
By default, this is implemented as a simple passthrough to the solver.
"""
(sfcc::SFCC)(soil::Soil, heat::Heat, state) = sfcc.solver(soil, heat, state, sfcc.f, sfcc.∇f)

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
Specialized implementation of Newton's method with backtracking line search for resolving
the energy conservation law, H = TC + Lθ. Attempts to find the root of the corresponding
temperature residual: ϵ = T - (H - Lθ(T)) / C(θ(T)) and uses backtracking to avoid
jumping over the solution. This prevents convergence issues that arise due to discontinuities
and non-monotonic behavior in most common soil freeze curves.
"""
@with_kw struct SFCCNewtonSolver <: SFCCSolver
    maxiter::Int = 10 # maximum number of iterations
    tol::float64 = 0.01 # absolute tolerance for convergence
    α₀::Float64 = 1.0 # initial step size multiplier
    τ::Float64 = 0.5 # step size decay for backtracking
end
function (s::SFCCNewtonSolver)(soil::Soil, heat::Heat{u"J"}, state, f, ∇f)
    # helper functions for handling arguments to freeze curve function, f;
    # some arguments may be grid state variables and thus we need to choose
    # the value at cell i; if not an array, we can just return x
    @inline atindex(x, i) = x
    @inline atindex(x::AbstractArray, i) = x[i]
    # helper function for updating θl, C, and the residual.
    @closure function residual(T, H, θw, θm, θo, f_args)
        θl = f(f_args...) # freeze curve
        C = heatcapacity(soil.hcparams, θw, θl, θm, θo)
        Tres = T - (H - θl*L) / C
        return Tres, θl, C
    end
    # get f arguments; note that this does create some redundancy in the arguments
    # eventually passed to the `residual` function; this is less than ideal but
    # probably shouldn't incur too much of a performance hit, just a few extra stack pointers!
    f_args = params(f, soil, heat, state)
    # iterate over each cell and solve the conservation law: H = TC + Lθ
    @inbounds @fastmath for i in 1:length(state.grids.T)
        itercount = 0
        let T = state.T[i], # temperature
            H = state.H[i], # enthalpy
            C = state.C[i], # heat capacity
            θl = state.θl[i], # liquid water content
            θtot = state.θw[i], # total water content
            θm = state.θm[i], # mineral content
            θo = state.θo[i], # organic content
            θp = state.θp[i], # porosity and/or θsat
            cw = soil.hcparams.cw, # heat capacity of liquid water
            f_argsᵢ = map(atindex, f_args);
            # compute initial residual
            Tres, θl, C = residual(T, H, θw, θm, θo, f_argsᵢ)
            while abs(Tres) > s.tol
                # derivative of freeze curve
                ∂θ∂T = ∇f(f_argsᵢ...)
                # derivative of residual by quotient rule;
                # note that this assumes heatcapacity to be a simple weighted average!
                # in the future, it might be a good idea to compute an automatic derivative
                # of heatcapacity in addition to the freeze curve function.
                ∂Tres∂T = 1.0 - ∂θ∂T*(-L*C - (H - θ*L)*cw)/C^2
                α = s.α₀ / ∂Tres∂T
                T̂ = T - α*Tres
                # do first residual check outside of loop;
                # this way, we don't decrease α unless we have to.
                T̂res, θl, C = residual(T̂, H, θw, θm, θo, f_argsᵢ)
                while sign(T̂res) != sign(Tres)
                    α = α*s.τ # decrease step size by τ
                    T̂ = T - α*Tres # new guess for T
                    T̂res, θl, C = residual(T̂, H, θw, θm, θo, f_argsᵢ)
                end
                T = T̂ # update T
                Tres = T̂res # update residual
                itercount += 1
                if itercount > maxiter
                    error("grid cell $i failed to converge after $maxiter iterations; residual: $(Tres[i])")
                end
            end
            # update state variables for cell i
            state.T[i] = T
            state.θl[i] = θl
            state.C[i] = C
        end
    end
    nothing
end
