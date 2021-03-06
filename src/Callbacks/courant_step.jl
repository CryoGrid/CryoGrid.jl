struct CFLHeatState{A}
    Δt::A
    default_dt::Float64
end
function (fn::CryoGridCallbackFunction{<:CFLHeatState{<:AbstractArray}})(u,p,t)
    let Δt = fn.state.Δt,
        Δx = dustrip(Δ(fn.setup.grid)),
        Ceff = getvar(:Ceff, fn.setup, u), # apparent heat capacity
        kc = getvar(:kc, fn.setup, u); # thermal cond. at grid centers
        @. Δt = Utils.adstrip(Δx^2 * Ceff / kc)
        Δt_min = minimum(Δt)
        IfElse.ifelse(isfinite(Δt_min) && Δt_min > 0, Δt_min, fn.state.default_dt)
    end
end
function (fn::CryoGridCallbackFunction{<:CFLHeatState{Nothing}})(u,p,t)
    let Δx = dustrip(Δ(fn.setup.grid)),
        Ceff = getvar(:Ceff, fn.setup, u), # apparent heat capacity
        kc = getvar(:kc, fn.setup, u); # thermal cond. at grid centers
        Δt = Utils.adstrip(Δx^2 * Ceff / kc)
        Δt_min = minimum(Δt)
        IfElse.ifelse(isfinite(Δt_min) && Δt_min > 0, Δt_min, fn.state.default_dt)
    end
end
function CFLStepLimiter(setup::HeatOnlySetup; courant_number=1//2, default_dt=60.0, iip=true)
    state = iip ? CFLHeatState(zero(dustrip(Δ(setup.grid))), default_dt) : CFLHeatState(nothing, default_dt)
    fn = CryoGridCallbackFunction(setup, state)
    StepsizeLimiter(fn; safety_factor=courant_number, max_step=true)
end

export CFLStepLimiter