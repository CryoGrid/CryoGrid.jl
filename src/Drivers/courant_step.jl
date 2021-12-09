struct CFLStepLimiter{TTile,A}
    tile::TTile
    Δt::A
    default_dt::Float64
end
function (cfl::CFLStepLimiter{<:Tile,<:AbstractArray})(u,p,t)
    let Δt = cfl.Δt,
        Δx = dustrip(Δ(cfl.tile.grid)),
        Ceff = getvar(:Ceff, cfl.tile, u), # apparent heat capacity
        kc = getvar(:kc, cfl.tile, u); # thermal cond. at grid centers
        @. Δt = Utils.adstrip(Δx^2 * Ceff / kc)
        Δt_min = minimum(Δt)
        IfElse.ifelse(isfinite(Δt_min) && Δt_min > 0, Δt_min, cfl.default_dt)
    end
end
function (cfl::CFLStepLimiter{<:Tile,Nothing})(u,p,t)
    let Δt = cfl.Δt,
        Δx = dustrip(Δ(cfl.tile.grid)),
        Ceff = getvar(:Ceff, cfl.tile, u), # apparent heat capacity
        kc = getvar(:kc, cfl.tile, u); # thermal cond. at grid centers
        Δt = Utils.adstrip(Δx^2 * Ceff / kc)
        Δt_min = minimum(Δt)
        IfElse.ifelse(isfinite(Δt_min) && Δt_min > 0, Δt_min, cfl.default_dt)
    end
end
function CFLStepLimiter(tile::HeatOnlyTile; courant_number=1//2, default_dt=60.0, iip=true)
    cfl = iip ? CFLStepLimiter(tile, zero(dustrip(Δ(setup.grid))), default_dt) : CFLHeatState(tile, nothing, default_dt)
    StepsizeLimiter(cfl; safety_factor=courant_number, max_step=true)
end
