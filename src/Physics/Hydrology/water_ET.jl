"""
    EvapTop <: Evapotranspiration

Represents a simple evaporation-only scheme where water is drawn only from the top-most grid cell.
Corresponds to evapotranspiration scheme 3 described in section 2.2.4 of Westermann et al. (2022).
"""
struct EvapTop <: Evapotranspiration end
"""
    DampedET{Tftr,Tfev,Tdtr,Tdev}

Corresponds to evapotranspiration scheme 2 described in section 2.2.4 of Westermann et al. (2022).
"""
Base.@kwdef struct DampedET{Tftr,Tdtr,Tdev} <: Evapotranspiration
    f_tr::Tftr = 0.5
    d_tr::Tdtr = 0.5u"m"
    d_ev::Tdev = 0.1u"m"
end
CryoGrid.parameterize(et::DampedET) = DampedET(
    f_tr = CryoGrid.parameterize(et.f_tr, domain=0..1, desc="Factor between 0 and 1 weighting transpirative vs. evaporative fluxes."),
    d_tr = CryoGrid.parameterize(et.d_tr, domain=0..Inf, desc="Damping depth for transpiration."),
    d_ev = CryoGrid.parameterize(et.d_ev, domain=0..Inf, desc="Damping depth for evaporation."),
)
"""
    evapotranspiration!(::SubSurface, ::WaterBalance, state)

Computes diagnostic evapotranspiration quantities for the given layer and water balance configuration, storing the results in `state`.
This method should generally be called *before* `interact!` for `WaterBalance`, e.g. in `diagnosticstep!`.
"""
function evapotranspiration!(::SubSurface, ::WaterBalance, state) end
function evapotranspiration!(
    sub::SubSurface,
    water::WaterBalance{<:BucketScheme,<:DampedET},
    state
)
    Δz = Δ(state.grid)
    z = cells(state.grid)
    et = water.et
    @inbounds for i in eachindex(z)
        state.w_ev[i] = Δz[i]*exp(-z[i] / et.d_ev)
        state.w_tr[i] = Δz[i]*exp(-z[i] / et.d_tr)
        let θwi = state.θwi[i],
            θfc = minwater(sub, water);
            state.αᶿ[i] = ifelse(θwi < θfc, 0.25(1-cos(π*θwi/θfc))^2, one(θwi))
        end
    end
    @inbounds @. state.f_et = et.f_tr*(state.αᶿ*state.w_tr) + (1-et.f_tr)*(state.αᶿ*state.w_ev)
end
"""
    ETflux(::SubSurface, water::WaterBalance{<:WaterFlow,<:Evapotranspiration}, state)

Computes the ET base flux as `Qe / (Lsg*ρw)` where `state.Qe` is typically provided as a boundary condition.
"""
function ETflux(::SubSurface, water::WaterBalance{<:WaterFlow,<:Evapotranspiration}, state)
    Qe = getscalar(state.Qe)
    Lsg = water.prop.Lsg # specific latent heat of vaporization
    ρw = water.prop.ρw
    return Qe / (Lsg*ρw)
end
"""
    evapotranspirative_fluxes!(::SubSurface, ::WaterBalance, state)

Computes diagnostic evapotranspiration quantities for the given layer and water balance configuration, storing the results in `state`.
This method should generally be called *after* `interact!` for `WaterBalance`, e.g. in `prognosticstep!`.
"""
function evapotranspirative_fluxes!(sub::SubSurface, water::WaterBalance, state) end
function evapotranspirative_fluxes!(
    sub::SubSurface,
    water::WaterBalance{<:BucketScheme,<:DampedET},
    state
)
    f_norm = sum(parent(state.f_et))
    # I guess we just ignore the flux at the lower boundary here... it will either be set
    # by the next layer or default to zero if no evapotranspiration occurs in the next layer.
    @inbounds for i in eachindex(cells(state.grid))
        fᵢ = IfElse.ifelse(f_norm > zero(f_norm), state.f_et[i] / f_norm, 0.0)
        state.jwET[i] += fᵢ * ETflux(sub, water, state)
        # add ET fluxes to total water flux
        state.jw[i] += state.jwET[i]
    end
end
# allow top evaporation-only scheme to apply by default for any water flow scheme
function evapotranspirative_fluxes!(sub::SubSurface, water::WaterBalance{<:WaterFlow,EvapTop}, state)
    state.jwET[1] += ETflux(sub, water, state)
    # add ET fluxes to total water flux
    state.jw[1] += state.jwET[1]
end
# CryoGrid methods
CryoGrid.basevariables(::Evapotranspiration) = (
    Diagnostic(:Qe, Scalar, u"J/s/m^2", desc="Latent heat flux at the surface."), # must be supplied by an interaction
)
CryoGrid.variables(et::DampedET) = (
    CryoGrid.basevariables(et)...,
    Diagnostic(:f_et, OnGrid(Cells), u"m", domain=0..1, desc="Evapotranspiration reduction factor."),
    Diagnostic(:w_ev, OnGrid(Cells), u"m", desc="Damped grid cell weight for evaporation."),
    Diagnostic(:w_tr, OnGrid(Cells), u"m", desc="Damped grid cell weight for transpiration"),
    Diagnostic(:αᶿ, OnGrid(Cells), domain=0..1, desc="Water availability coefficient."),
)
function CryoGrid.interact!(
    ::SubSurface,
    ::WaterBalance{<:BucketScheme,<:DampedET},
    ::SubSurface, ::WaterBalance{<:BucketScheme,<:DampedET},
    state1,
    state2
)
    # propagate surface latent heat flux to next layer
    state2.Qe .= state1.Qe
end
