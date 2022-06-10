using NLsolve

"""
    SFCCNewtonSolver <: SFCCSolver

Specialized implementation of Newton's method with backtracking line search for resolving
the energy conservation law, H = TC + Lθ. Attempts to find the root of the corresponding
temperature residual: ϵ = T - (H - Lθ(T)) / C(θ(T)) and uses backtracking to avoid
jumping over the solution. This prevents convergence issues that arise due to
discontinuities and strong non-linearity in most common soil freeze curves.
"""
Base.@kwdef struct SFCCNewtonSolver <: SFCCSolver
    maxiter::Int = 100 # maximum number of iterations
    abstol::Float64 = 1e-2 # absolute tolerance for convergence
    reltol::Float64 = 1e-4 # relative tolerance for convergence
    α₀::Float64 = 1.0 # initial step size multiplier
    τ::Float64 = 0.7 # step size decay for backtracking
end
# Helper function for updating θw, C, and the residual.
@inline function sfccresidual(soil::Soil, heat::Heat, f::F, f_args::Fargs, f_hc, T, H) where {F,Fargs}
    L = heat.prop.L
    θw = f(T, f_args...)
    C = f_hc(θw)
    Tres = T - (H - θw*L) / C
    return Tres, θw, C
end
# Newton solver implementation
function sfccsolve(solver::SFCCNewtonSolver, soil::Soil, heat::Heat, f, f_args, f_hc, H, T₀::Nothing=nothing)
    T₀ = IfElse.ifelse(H < zero(H), H / f_hc(0.0), zero(H))
    return sfccsolve(solver, soil, heat, f, f_args, f_hc, H, T₀)
end
using ForwardDiff
function sfccsolve(solver::SFCCNewtonSolver, soil::Soil, heat::Heat, f::F, f_args, f_hc, H, T₀) where {F}
    T = T₀
    L = heat.prop.L
    cw = heat.prop.cw # heat capacity of liquid water
    ci = heat.prop.ci # heat capacity of ice
    α₀ = solver.α₀
    τ = solver.τ
    # compute initial residual
    Tres, θw, C = sfccresidual(soil, heat, f, f_args, f_hc, T, H)
    itercount = 0
    T_converged = false
    while abs(Tres) > solver.abstol && abs(Tres) / abs(T) > solver.reltol
        if itercount >= solver.maxiter
            return (;T, Tres, θw, itercount, T_converged)
        end
        # derivative of freeze curve
        ∂θ∂T = ForwardDiff.derivative(T -> f(T, f_args...), T)
        # derivative of residual by quotient rule;
        # note that this assumes heatcapacity to be a simple weighted average!
        # in the future, it might be a good idea to compute an automatic derivative
        # of heatcapacity in addition to the freeze curve function.
        ∂Tres∂T = 1.0 - ∂θ∂T*(-L*C - (H - θw*L)*(cw-ci))/C^2
        α = α₀ / ∂Tres∂T
        T̂ = T - α*Tres
        # do first residual check outside of loop;
        # this way, we don't decrease α unless we have to.
        T̂res, θw, C = sfccresidual(soil, heat, f, f_args, f_hc, T̂, H)
        inneritercount = 0
        # simple backtracking line search to avoid jumping over the solution
        while sign(T̂res) != sign(Tres)
            if inneritercount > 100
                @warn "Backtracking failed; this should not happen. Current state: α=$α, T=$T, T̂=$T̂, residual $(T̂res), initial residual: $(Tres)"
                return (;T, Tres, θw, itercount, T_converged)
            end
            α = α*τ # decrease step size by τ
            T̂ = T - α*Tres # new guess for T
            T̂res, θw, C = sfccresidual(soil, heat, f, f_args, f_hc, T̂, H)
            inneritercount += 1
        end
        T = T̂ # update T
        Tres = T̂res # update residual
        itercount += 1
    end
    T_converged = true
    return (;T, Tres, θw, itercount, T_converged)
end
function HeatConduction.enthalpyinv(soil::Soil, heat::Heat{<:SFCC{F,SFCCNewtonSolver},Enthalpy}, state, i) where {F}
    sfcc = freezecurve(heat)
    f = sfcc.f
    # get f arguments; note that this does create some redundancy in the arguments
    # eventually passed to the `residual` function; this is less than ideal but
    # probably shouldn't incur too much of a performance hit, just a few extra stack pointers!
    f_args = sfccargs(f, soil, heat, state)
    solver = sfcc.solver
    @inbounds let T₀ = i > 1 ? state.T[i-1] : nothing,
        H = state.H[i] |> Utils.adstrip, # enthalpy
        f_hc = partial(heatcapacity, liquidwater, soil, heat, state, i),
        f_argsᵢ = Utils.selectat(i, Utils.adstrip, f_args);
        T, _, _, _, _ = sfccsolve(solver, soil, heat, f, f_argsᵢ, f_hc, H, T₀)
        return T
    end
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
    errtol::Float64 | false
    SFCCPreSolver(cache, Tmin, errtol) = new{typeof(cache)}(cache, Tmin, errtol)
    """
        SFCCPreSolver(Tmin=-60.0, errtol=1e-4)

    Constructs a new `SFCCPreSolver` with minimum temperature `Tmin` and integration step `dH`.
    Enthalpy values below `H(Tmin)` under the given freeze curve will be extrapolated with a
    constant/flat function. `errtol` determines the permitted local error in the interpolant.
    """
    function SFCCPreSolver(;Tmin=-60.0, errtol=1e-4)
        cache = SFCCPreSolverCache()
        new{typeof(cache)}(cache, Tmin, errtol)
    end
end
mutable struct SFCCPreSolverCache{F,∇F}
    f::F # H⁻¹ interpolant
    ∇f::∇F # derivative of f
    function SFCCPreSolverCache()
        # initialize with dummy functions to get type information
        x  = -3e8:1e6:3e8
        dummy_f = _build_interpolant(x, zeros(length(x)))
        dummy_∇f = first ∘ _interpgrad(dummy_f)
        return new{typeof(dummy_f),typeof(dummy_∇f)}(dummy_f, dummy_∇f)
    end
end
_interpgrad(f) = (args...) -> Interpolations.gradient(f, args...)
function _build_interpolant(Hs, θs)
    return Interpolations.extrapolate(
        Interpolations.interpolate((Vector(Hs),), θs, Interpolations.Gridded(Interpolations.Linear())),
        Interpolations.Flat()
    )
end
function CryoGrid.initialcondition!(soil::Soil{<:HomogeneousCharacteristicFractions}, heat::Heat, sfcc::SFCC{F,<:SFCCPreSolver}, state) where {F}
    L = heat.prop.L
    args = sfccargs(sfcc.f, soil, heat, state)
    state.θw .= sfcc.f.(state.T, args...)
    heatcapacity!(soil, heat, state)
    @. state.H = enthalpy(state.T, state.C, L, state.θw)
    # pre-solve freeze curve;
    # note that this is only valid given that the following assumptions hold:
    # 1) none of the freeze curve parameters (e.g. soil properties) change
    # 2) soil properties are uniform in `soil`
    let θsat = args[1],
        θtot = args[2],
        tail_args = args[3:end],
        θm = mineral(soil, state),
        θo = organic(soil, state),
        θa = 1-θm-θo-θtot,
        L = heat.prop.L,
        Tmin = sfcc.solver.Tmin,
        Tmax = 0.0,
        f(T) = sfcc.f(T, θsat, θtot, tail_args...),
        C(θw) = heatcapacity(soil, heat, θw, θtot-θw, θa, θm, θo),
        Hmin = enthalpy(Tmin, C(f(Tmin)), L, f(Tmin)),
        Hmax = enthalpy(Tmax, C(f(Tmax)), L, f(Tmax));
        sfccsolver = SFCCNewtonSolver(abstol=sfcc.solver.errtol, reltol=sfcc.solver.errtol, maxiter=100)
        # residual as a function of T and H
        resid(T,H) = sfccresidual(soil, heat, sfcc.f, args, C, T, H)
        function solve(H,T₀)
            opt = sfccsolve(sfccsolver, soil, heat, sfcc.f, args, C, H, T₀)
            @assert opt.T_converged "solver failed to converge after $(opt.itercount) iterations: H=$H, T₀=$T₀, T=$(opt.T), Tres=$(opt.Tres), θw=$(opt.θw)"
            return opt
        end
        function deriv(T) # implicit partial derivative w.r.t H as a function of T
            θw, ∂θ∂T = ∇(f, T)
            # get C_eff, i.e. dHdT
            ∂H∂T = HeatConduction.C_eff(T, C(θw), heat.prop.L, ∂θ∂T, heat.prop.cw, heat.prop.ci)
            # valid by chain rule and inverse function theorem
            return ∂θ∂T / ∂H∂T
        end
        function step(ΔH, H, θ, ∂θ∂H, T₀)
            # get first order estimate
            θest = θ + ΔH*∂θ∂H
            # get true θ at H + ΔH
            θsol = solve(H + ΔH, T₀).θw
            err = abs(θsol - θest)
            # return residual of error with target error
            return err
        end
        T = [Tmin]
        H = [Hmin]
        θ = [f(T[1])]
        ∂θ∂H = [deriv(T[1])]
        @assert isfinite(H[1]) && isfinite(θ[1]) "H=$H, θ=$θ"
        while H[end] < Hmax
            # find the optimal step size
            ϵ = Inf
            ΔH = heat.prop.L*θtot/10 # initially set to large value
            while abs(ϵ) > sfcc.solver.errtol
                ϵ = step(ΔH, H[end], θ[end], ∂θ∂H[end], T[end])
                # iteratively halve the step size until error tolerance is satisfied
                ΔH *= 0.5
            end
            Hnew = H[end] + ΔH
            @assert isfinite(Hnew) "isfinite(ΔH) failed; H=$(H[end]), T=$(T[end]), ΔH=$ΔH"
            opt = solve(Hnew, T[end])
            push!(H, Hnew)
            push!(θ, opt.θw)
            push!(T, opt.T)
            push!(∂θ∂H, deriv(opt.T))
        end
        sfcc.solver.cache.f = _build_interpolant(H, θ)
        sfcc.solver.cache.∇f = first ∘ _interpgrad(sfcc.solver.cache.f)
    end
end
