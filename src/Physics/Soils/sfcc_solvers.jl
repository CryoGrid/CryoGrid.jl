"""
    SFCCNewtonSolver <: SFCCSolver

Specialized implementation of Newton's method with backtracking line search for resolving
the energy conservation law, H = TC + Lθ. Attempts to find the root of the corresponding
temperature residual: ϵ = T - (H - Lθ(T)) / C(θ(T)) and uses backtracking to avoid
jumping over the solution. This prevents convergence issues that arise due to
discontinuities and strong non-linearity in most common soil freeze curves.
"""
@with_kw struct SFCCNewtonSolver <: SFCCSolver
    maxiter::Int = 100 # maximum number of iterations
    abstol::Float64 = 1e-2 # absolute tolerance for convergence
    reltol::Float64 = 1e-2 # relative tolerance for convergence
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
function sfccsolve(solver::SFCCNewtonSolver, soil::Soil, f, ∇f, f_args, H, L, θw, θm, θo, T₀=nothing)
    # compute initial guess T by setting θl according to free water scheme
    T = if isnothing(T₀)
        let Lθ = L*θw;
            if H < 0
                H / heatcapacity(soil,θw,0.0,θm,θo)
            elseif H >= 0 && H < Lθ
                (1.0 - H/Lθ)*0.1
            else
                (H - Lθ) / heatcapacity(soil,θw,θw,θm,θo)
            end
        end
    else
        T₀
    end
    cw = soil.hc.cw # heat capacity of liquid water
    α₀ = solver.α₀
    τ = solver.τ
    # compute initial residual
    Tres, θl, C = residual(T, H, θw, θm, θo, L, soil, f, f_args)
    itercount = 0
    @fastmath while abs(Tres) > solver.abstol && abs(Tres) / abs(T) > solver.reltol
        if itercount > solver.maxiter
            convergencefailure(solver.onfail, i, solver.maxiter, Tres)
            break
        end
        # derivative of freeze curve
        args = tuplejoin((T,),f_args)
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
        T̂res, θl, C = residual(T̂, H, θw, θm, θo, L, soil, f, f_args)
        inneritercount = 0
        # simple backtracking line search to avoid jumping over the solution
        while sign(T̂res) != sign(Tres)
            if inneritercount > 100
                @warn "Backtracking failed; this should not happen. Current state: α=$α, T=$T, T̂=$T̂, residual $(T̂res), initial residual: $(Tres)"
                break
            end
            α = α*τ # decrease step size by τ
            T̂ = T - α*Tres # new guess for T
            T̂res, θl, C = residual(T̂, H, θw, θm, θo, L, soil, f, f_args)
            inneritercount += 1
        end
        T = T̂ # update T
        Tres = T̂res # update residual
        itercount += 1
    end
    return (;T, Tres, θl, itercount)
end
function (solver::SFCCNewtonSolver)(soil::Soil, heat::Heat{<:SFCC,Enthalpy}, state, f, ∇f)
    # get f arguments; note that this does create some redundancy in the arguments
    # eventually passed to the `residual` function; this is less than ideal but
    # probably shouldn't incur too much of a performance hit, just a few extra stack pointers!
    f_args = sfccparams(f, soil, heat, state)
    # iterate over each cell and solve the conservation law: H = TC + Lθ
    @threaded for i in 1:length(state.T)
        @inbounds @fastmath let T₀ = i > 1 ? state.T[i-1] : nothing,
            H = state.H[i] |> Utils.adstrip, # enthalpy
            θl = state.θl[i] |> Utils.adstrip, # liquid water content
            θw = totalwater(soil, heat, state, i) |> Utils.adstrip, # total water content
            θm = mineral(soil, heat, state, i) |> Utils.adstrip, # mineral content
            θo = organic(soil, heat, state, i) |> Utils.adstrip, # organic content
            L = heat.L, # specific latent heat of fusion
            f_argsᵢ = Utils.selectat(i, Utils.adstrip, f_args);
            T, _, _, _ = sfccsolve(solver, soil, f, ∇f, f_argsᵢ, H, L, θw, θm, θo, T₀)
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
                    state.dHdT[i] = state.C[i] + dθdT
                    state.T[i] = (H - L*θl) / state.C[i]
                end
            end
        end
    end
    return nothing
end
"""
    SFCCPreSolver{TCache} <: SFCCSolver

A fast SFCC "solver" which pre-builds an interpolant for the freeze curve in terms of enthalpy, θ(H).
Note that this solver is **only valid when all freeze curve parameters are held constant** and will
produce incorrect results otherwise.
"""
@flattenable struct SFCCPreSolver{TCache} <: SFCCSolver
    cache::TCache | false
    Tmin::Float64 | false
    dH::Float64 | false
    SFCCPreSolver(cache, Tmin, dH) = new{typeof(cache)}(cache, Tmin, dH)
    """
        SFCCPreSolver(Tmin=-50.0, dH=2e5)

    Constructs a new `SFCCPreSolver` with minimum temperature `Tmin` and integration step `dH`.
    Enthalpy values below `H(Tmin)` under the given freeze curve will be extrapolated with a
    constant/flat function. `dH` determines the step size used when integrating `dθdH`; smaller
    values will produce a more accurate interpolant at the cost of storing more knots and slower
    initialization. The default value of `dH=2e5` should be sufficient for most use-cases.
    """
    function SFCCPreSolver(;Tmin=-50.0, dH=2e5)
        cache = SFCCPreSolverCache()
        new{typeof(cache)}(cache, Tmin, dH)
    end
end
mutable struct SFCCPreSolverCache
    f # H⁻¹ interpolant
    SFCCPreSolverCache() = new()
end
function initialcondition!(soil::Soil{<:HomogeneousCharacteristicFractions}, heat::Heat, sfcc::SFCC{F,∇F,<:SFCCPreSolver}, state) where {F,∇F}
    L = heat.L
    params = sfccparams(sfcc.f, soil, heat, state)
    state.θl .= sfcc.f.(state.T, params...)
    heatcapacity!(soil, heat, state)
    @. state.H = enthalpy(state.T, state.C, L, state.θl)
    # pre-solve freeze curve;
    # note that this is only valid given that the following assumptions hold:
    # 1) none of the freeze curve parameters (e.g. soil properties) change
    # 2) soil properties are uniform in `soil`
    let Tₘ = params[1],
        θres = params[2],
        θsat = params[3],
        θtot = params[4],
        args = params[5:end],
        θm = mineral(soil, state),
        θo = organic(soil, state),
        L = heat.L,
        Tmin = sfcc.solver.Tmin,
        Tmax = Tₘ,
        θ(T) = sfcc.f(T, Tₘ, θres, θsat, θtot, args...),
        C(T) = heatcapacity(soil, θtot, θ(T), θm, θo),
        Hmin = enthalpy(Tmin, C(Tmin), L, θ(Tmin)),
        Hmax = enthalpy(Tmax, C(Tmax), L, θ(Tmax)),
        dH = sfcc.solver.dH,
        Hs = Hmin:dH:Hmax;
        θs = Vector{eltype(state.θl)}(undef, length(Hs))
        θs[1] = θ(Tmin)
        Ts = Vector{eltype(state.T)}(undef, length(Hs))
        Ts[1] = Tmin
        solver = SFCCNewtonSolver(abstol=1e-3, reltol=1e-3, onfail=:error)
        for i in 2:length(Hs)
            Hᵢ = Hs[i]
            T₀ = Ts[i-1] # use previous temperature value as initial guess
            res = sfccsolve(solver, soil, sfcc.f, sfcc.∇f, params, Hᵢ, L, θtot, θm, θo, T₀)
            θs[i] = res.θl
            Ts[i] = res.T
        end
        sfcc.solver.cache.f = Interpolations.extrapolate(
            Interpolations.interpolate((Vector(Hs),), θs, Interpolations.Gridded(Interpolations.Linear())),
            Interpolations.Flat()
        )
    end
end
function (s::SFCCPreSolver)(soil::Soil{<:HomogeneousCharacteristicFractions}, heat::Heat, state, _, _)
    state.θl .= s.cache.f.(state.H)
    heatcapacity!(soil, heat, state)
    @. state.T = (state.H - heat.L*state.θl) / state.C
    ∇f = first ∘ ∇(s.cache.f)
    @. state.dHdT = 1 / ∇f.(state.H)
end
