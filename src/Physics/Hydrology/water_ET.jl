"""
    DampedET{Tftr,Tfev,Tdtr,Tdev}

Corresponds to evapotranspiration scheme 3 described in section 2.2.2 of Westermann et al. (2022).
"""
Base.@kwdef struct DampedET{Tftr,Tdtr,Tdev} <: Evapotranspiration
    f_tr::Tftr = Param(0.5)
    d_tr::Tdtr = Param(0.5, units=u"m")
    d_ev::Tdev = Param(0.1, units=u"m")
end
"""
    evapotranspiration!(::SubSurface, ::WaterBalance, state)

Computes evapotranspiration fluxes for the given layer and water balance configuration, storing the results in `state`.
This method should generally be called after `interact!` for `WaterBalance`, e.g. in `prognosticstep!`.
"""
evapotranspiration!(::SubSurface, ::WaterBalance, state) = nothing
function evapotranspiration!(sub::SubSurface, water::WaterBalance{<:BucketScheme,<:DampedET}, state)
    Δz = Δ(state.grid)
    z = cells(state.grid)
    et = water.et
    @inbounds for i in eachindex(Δz)
        state.w_ev[i] = Δz[i]*exp(-z[i] / et.d_ev)
        state.w_tr[i] = Δz[i]*exp(-z[i] / et.d_tr)
        let θwi = state.θwi[i],
            θfc = fieldcapacity(sub, water);
            state.αᶿ[i] = ifelse(θwi < θfc, 0.25(1-cos(π*θwi/θfc))^2, one(θwi))
        end
    end
    @inbounds @. state.f_et = et.f_tr*(state.αᶿ*state.w_tr) + (1-et.f_tr)*(state.αᶿ*state.w_ev)
    norm_f = sum(state.f_et)
    Lsg = water.prop.consts.Lsg # specific latent heat of vaporization
    ρw = water.prop.consts.ρw
    # I guess we just ignore the flux at the lower boundary here... it will either be set
    # by the next layer or default to zero if no evapotranspiration occurs in the next layer.
    Qe = getscalar(state.Qe)
    @inbounds for i in eachindex(Δz)
        state.jwET[i] += state.f_et[i] / norm_f * Qe / (Lsg*ρw)
        # add ET fluxes to total water flux
        state.jw[i] += state.jwET[i]
    end
end
CryoGrid.basevariables(::Evapotranspiration) = (
    Diagnostic(:Qe, Scalar, desc="Latent heat flux at the surface."), # must be supplied by an interaction
)
CryoGrid.variables(et::DampedET) = (
    CryoGrid.basevariables(et)...,
    Diagnostic(:f_et, OnGrid(Cells), domain=0..1, desc="Evapotranspiration reduction factor."),
    Diagnostic(:w_ev, OnGrid(Cells), desc="Damped grid cell weight for evaporation."),
    Diagnostic(:w_tr, OnGrid(Cells), desc="Damped grid cell weight for transpiration"),
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
