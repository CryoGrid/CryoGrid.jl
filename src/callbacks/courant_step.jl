struct CFLHeatState{A}
    Δt::A
end
function (fn::CryoGridCallbackFunction{<:CFLHeatState})(u,p,t)
    let Δt = fn.state.Δt,
        Δx = Δ(fn.setup.grid),
        Ceff = getvar(:Ceff, fn.setup, u), # apparent heat capacity
        kc = getvar(:kc, fn.setup, u); # thermal cond. at grid centers
        @. Δt = Δx^2 * Ceff / kc
        minimum(Δt)
    end
end
function CFLStepLimiter(setup::HeatOnlySetup; courant_number=1//2)
    state = CFLHeatState(similar(Δ(setup.grid)))
    fn = CryoGridCallbackFunction(setup, state)
    StepsizeLimiter(fn; safety_factor=courant_number, max_step=true)
end

export CFLStepLimiter