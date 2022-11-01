# Thermal properties (generic)
"""
    enthalpy(T, C, L, θ) = T*C + L*θ

Discrete enthalpy function on temperature, heat capacity, specific latent heat of fusion, and liquid water content.
"""
@inline enthalpy(T, C, L, θ) = T*C + L*θ
"""
    enthalpyinv(H, C, L, θ) = (H - L*θ) / C

Discrete inverse enthalpy function given H, C, L, and θ.
"""
@inline enthalpyinv(H, C, L, θ) = (H - L*θ) / C
"""
    C_eff(T, C, L, ∂θw∂T, cw, ci) = C + ∂θw∂T*(L + T*(cw - ci))

Computes the apparent or "effective" heat capacity `∂H∂T` as a function of temperature, volumetric heat capacity,
latent heat of fusion, derivative of the freeze curve `∂θw∂T`, and the constituent heat capacities of water and ice.
"""
@inline C_eff(T, C, L, ∂θw∂T, cw, ci) = C + ∂θw∂T*(L + T*(cw - ci))
"""
    thermalconductivities(::SubSurface, heat::HeatBalance)

Get thermal conductivities for generic `SubSurface` layer.
"""
@inline function thermalconductivities(sub::SubSurface)
    @unpack kw, ki, ka = thermalproperties(sub)
    return (kw, ki, ka)
end
"""
    heatcapacities(::SubSurface, heat::HeatBalance)

Get heat capacities for generic `SubSurface` layer.
"""
@inline function heatcapacities(sub::SubSurface)
    @unpack cw, ci, ca = thermalproperties(sub)
    return cw, ci, ca
end

# Generic heat conduction implementation
"""
    freezethaw!(sub::SubSurface, heat::HeatBalance, state)

Calculates freezing and thawing effects, including evaluation of the freeze curve.
In general, this function should compute at least the liquid/frozen water contents
and the corresponding heat capacity. Other variables such as temperature or enthalpy
should also be computed depending on the thermal scheme being implemented.
"""
freezethaw!(::SubSurface, ::HeatBalance, state) = error("missing implementation of freezethaw!")
"""
    heatcapacity(sub::SubSurface, θfracs...)

Computes the heat capacity as a weighted average over constituent capacities with volumetric fractions `θfracs`.
"""
@inline function heatcapacity(sub::SubSurface, θfracs...)
    cs = heatcapacities(sub)
    return sum(map(*, cs, θfracs))
end
"""
    heatcapacity!(sub::SubSurface, state)

Computes the heat capacity for the given layer from the current state and stores the result in-place in the state variable `C`.
"""
@inline function heatcapacity!(sub::SubSurface, state)
    @inbounds for i in 1:length(state.T)
        θfracs = volumetricfractions(sub, state, i)
        state.C[i] = heatcapacity(sub, θfracs...)
    end
end
"""
    thermalconductivity(sub::SubSurface, θfracs...)

Computes the thermal conductivity as a squared weighted sum over constituent conductivities with volumetric fractions `θfracs`.
"""
@inline function thermalconductivity(sub::SubSurface, θfracs...)
    ks = thermalconductivities(sub)
    return sum(map(*, map(sqrt, ks), θfracs))^2
end
"""
    thermalconductivity!(sub::SubSurface, state)

Computes the thermal conductivity for the given layer from the current state and stores the result in-place in the state variable `C`.
"""
@inline function thermalconductivity!(sub::SubSurface, state)
    @inbounds for i in 1:length(state.T)
        θfracs = volumetricfractions(sub, state, i)
        state.kc[i] = thermalconductivity(sub, θfracs...)
    end
    # thermal conductivity at boundaries
    # assumes boundary conductivities = cell conductivities
    @inbounds state.k[1] = state.kc[1]
    @inbounds state.k[end] = state.kc[end]
    # Harmonic mean of inner conductivities
    @inbounds let k = (@view state.k[2:end-1]),
        Δk = Δ(state.grids.k);
        harmonicmean!(k, state.kc, Δk)
    end
end
"""
Variable definitions for heat conduction (enthalpy) on any SubSurface layer.
"""
CryoGrid.variables(heat::HeatBalance{<:FreezeCurve,<:Enthalpy}) = (
    Prognostic(:H, OnGrid(Cells), u"J/m^3"),
    Diagnostic(:T, OnGrid(Cells), u"°C"),
    CryoGrid.basevariables(heat)...,
)
"""
Variable definitions for heat conduction (temperature).
"""
CryoGrid.variables(heat::HeatBalance{<:FreezeCurve,<:Temperature}) = (
    Prognostic(:T, OnGrid(Cells), u"°C"),
    Diagnostic(:H, OnGrid(Cells), u"J/m^3"),
    Diagnostic(:∂H∂t, OnGrid(Cells), u"W/m^3"),
    CryoGrid.basevariables(heat)...,
)
"""
Common variable definitions for all heat implementations.
"""
CryoGrid.basevariables(::HeatBalance) = (
    Diagnostic(:jH, OnGrid(Edges), u"W/m^2"),    
    Diagnostic(:∂H∂T, OnGrid(Cells), u"J/K/m^3", domain=0..Inf),
    Diagnostic(:∂θw∂T, OnGrid(Cells), domain=0..Inf),
    Diagnostic(:C, OnGrid(Cells), u"J/K/m^3"),
    Diagnostic(:k, OnGrid(Edges), u"W/m/K"),
    Diagnostic(:kc, OnGrid(Cells), u"W/m/K"),
    Diagnostic(:θw, OnGrid(Cells), domain=0..1),
)
function resetfluxes!(sub::SubSurface, heat::HeatBalance, state)
    # Reset energy fluxes to zero; this is redundant when H is the prognostic variable
    # but necessary when it is not.
    @. state.∂H∂t = zero(eltype(state.∂H∂t))
    @. state.jH = zero(eltype(state.jH))
end
"""
Diagonstic step for heat conduction (all state configurations) on any subsurface layer.
"""
function CryoGrid.diagnosticstep!(sub::SubSurface, heat::HeatBalance, state)
    resetfluxes!(sub, heat, state)
    # Evaluate freeze/thaw processes
    freezethaw!(sub, heat, state)
    # Update thermal conductivity
    thermalconductivity!(sub, state)
    return nothing
end
# Boundary fluxes
@inline CryoGrid.boundaryflux(::Neumann, bc::HeatBC, top::Top, heat::HeatBalance, sub::SubSurface, stop, ssub) = boundaryvalue(bc,top,heat,sub,stop,ssub)
@inline CryoGrid.boundaryflux(::Neumann, bc::HeatBC, bot::Bottom, heat::HeatBalance, sub::SubSurface, sbot, ssub) = boundaryvalue(bc,bot,heat,sub,sbot,ssub)
@inline function CryoGrid.boundaryflux(::Dirichlet, bc::HeatBC, top::Top, heat::HeatBalance, sub::SubSurface, stop, ssub)
    Δk = CryoGrid.thickness(sub, ssub, first) # using `thickness` allows for generic layer implementations
    @inbounds let Tupper=boundaryvalue(bc,top,heat,sub,stop,ssub),
        Tsub=ssub.T[1],
        k=ssub.k[1],
        δ=Δk/2; # distance to boundary
        -k*(Tsub-Tupper)/δ
    end
end
@inline function CryoGrid.boundaryflux(::Dirichlet, bc::HeatBC, bot::Bottom, heat::HeatBalance, sub::SubSurface, sbot, ssub)
    Δk = CryoGrid.thickness(sub, ssub, last) # using `thickness` allows for generic layer implementations
    @inbounds let Tlower=boundaryvalue(bc,bot,heat,sub,sbot,ssub),
        Tsub=ssub.T[end],
        k=ssub.k[end],
        δ=Δk/2; # distance to boundary
        # note again the inverted sign; positive here means *upward from* the bottom boundary
        k*(Tlower-Tsub)/δ
    end
end
"""
Generic top interaction. Computes flux jH at top cell.
"""
function CryoGrid.interact!(top::Top, bc::HeatBC, sub::SubSurface, heat::HeatBalance, stop, ssub)
    # boundary flux
    ssub.jH[1] += boundaryflux(bc, top, heat, sub, stop, ssub)
    return nothing
end
"""
Generic bottom interaction. Computes flux jH at bottom cell.
"""
function CryoGrid.interact!(sub::SubSurface, heat::HeatBalance, bot::Bottom, bc::HeatBC, ssub, sbot)
    # boundary flux; here we flip the sign since a positive flux is by convention downward
    ssub.jH[end] -= boundaryflux(bc, bot, heat, sub, sbot, ssub)
    return nothing
end
"""
Generic subsurface interaction. Computes flux jH at boundary between subsurface layers.
"""
function CryoGrid.interact!(sub1::SubSurface, ::HeatBalance, sub2::SubSurface, ::HeatBalance, s1, s2)
    Δk₁ = CryoGrid.thickness(sub1, s1, last)
    Δk₂ = CryoGrid.thickness(sub2, s2, first)
    # thermal conductivity between cells
    k = s1.k[end] = s2.k[1] =
        @inbounds let k₁ = s1.kc[end],
            k₂ = s2.kc[1],
            Δ₁ = Δk₁[end],
            Δ₂ = Δk₂[1];
            harmonicmean(k₁, k₂, Δ₁, Δ₂)
        end
    # calculate heat flux between cells (positive downward)
    Qᵢ = @inbounds let z₁ = CryoGrid.midpoint(sub1, s1, last),
        z₂ = CryoGrid.midpoint(sub2, s2, first),
        δ = z₂ - z₁;
        -k*(s2.T[1] - s1.T[end]) / δ
    end
    # set inner boundary flux
    s1.jH[end] += Qᵢ
    # these should already be equal since jH[1] on s2 and jH[end] on s1 should point to the same array index;
    # but we set them explicitly just to be safe.
    s2.jH[1] = s1.jH[end]
    return nothing
end
"""
Prognostic step for heat conduction (enthalpy) on subsurface layer.
"""
function CryoGrid.prognosticstep!(::SubSurface, ::HeatBalance{<:FreezeCurve,<:Enthalpy}, state)
    Δk = Δ(state.grids.k) # cell sizes
    ΔT = Δ(state.grids.T) # midpoint distances
    # compute internal fluxes and non-linear diffusion assuming boundary fluxes have been set
    nonlineardiffusion!(state.∂H∂t, state.jH, state.T, ΔT, state.k, Δk)
    return nothing
end
"""
Prognostic step for heat conduction (temperature) on subsurface layer.
"""
function CryoGrid.prognosticstep!(sub::SubSurface, ::HeatBalance{<:FreezeCurve,<:Temperature}, state)
    Δk = Δ(state.grids.k) # cell sizes
    ΔT = Δ(state.grids.T) # midpoint distances
    # compute internal fluxes and non-linear diffusion assuming boundary fluxes have been set
    nonlineardiffusion!(state.∂H∂t, state.jH, state.T, ΔT, state.k, Δk)
    # Compute temperature flux by dividing by ∂H∂T;
    # ∂H∂T should be computed by the freeze curve.
    @inbounds @. state.∂T∂t = state.∂H∂t / state.∂H∂T
    return nothing
end
# CFL not defined for free-water freeze curve
CryoGrid.timestep(::SubSurface, heat::HeatBalance{FreeWater,<:Enthalpy,<:Physics.CFL}, state) = Inf
"""
    timestep(::SubSurface, ::HeatBalance{Tfc,THeatOp,CFL}, state) where {Tfc,THeatOp}

Implementation of `timestep` for `HeatBalance` using the Courant-Fredrichs-Lewy condition
defined as: Δt_max = u*Δx^2, where`u` is the "characteristic velocity" which here
is taken to be the diffusivity: `∂H∂T / kc`.
"""
function CryoGrid.timestep(::SubSurface, heat::HeatBalance{Tfc,THeatOp,<:Physics.CFL}, state) where {Tfc,THeatOp}
    derivative(::Enthalpy, state) = state.∂H∂t
    derivative(::Temperature, state) = state.∂T∂t
    prognostic(::Enthalpy, state) = state.H
    prognostic(::Temperature, state) = state.T
    Δx = Δ(state.grid)
    dtmax = Inf
    @inbounds for i in eachindex(Δx)
        dtmax = let v = state.∂H∂T[i] / state.kc[i], # characteristic velocity
            c = heat.dtlim.courant_number, # couant number
            Δx = Δx[i],
            Δt = c*v*Δx^2;
            # compute maxdelta timestep limit and choose the smaller one
            Δt = min(Δt, heat.dtlim.maxdelta(derivative(heat.op, state)[i], prognostic(heat.op, state)[i], state.t))
            # select smaller of current dtmax and Δt
            min(dtmax, Δt)
        end
    end
    dtmax = isfinite(dtmax) ? dtmax : Inf
    return dtmax
end
function CryoGrid.timestep(::SubSurface, heat::HeatBalance{Tfc,THeatOp,<:Physics.MaxDelta}, state) where {Tfc,THeatOp}
    Δx = Δ(state.grid)
    dtmax = Inf
    @inbounds for i in eachindex(Δx)
        dtmax = min(dtmax, heat.dtlim(state.∂H∂t[i], state.H[i], state.t))
    end
    dtmax = isfinite(dtmax) && dtmax > 0 ? dtmax : Inf
    return dtmax
end
# Free water freeze curve
@inline function enthalpyinv(sub::SubSurface, heat::HeatBalance{FreeWater,<:Enthalpy}, state, i)
    hc = partial(heatcapacity, Val{:θw}(), sub, state, i)
    return enthalpyinv(heat.freezecurve, hc, state.H[i], hc.θwi, heat.prop.L)
end
@inline function enthalpyinv(::FreeWater, hc::F, H, θwi, L) where {F}
    θw, I_t, I_f, I_c, Lθ = FreezeCurves.freewater(H, θwi, L)
    C = hc(θw)
    T = (I_t*(H-Lθ) + I_f*H)/C
    return T, θw, C
end
"""
    freezethaw!(sub::SubSurface, heat::HeatBalance{FreeWater,<:Enthalpy}, state)

Implementation of "free water" freezing characteristic for any subsurface layer.
Assumes that `state` contains at least temperature (T), enthalpy (H), heat capacity (C),
total water content (θwi), and liquid water content (θw).
"""
@inline function freezethaw!(sub::SubSurface, heat::HeatBalance{FreeWater,<:Enthalpy}, state)
    @inbounds for i in 1:length(state.H)
        # update T, θw, C
        state.T[i], state.θw[i], state.C[i] = enthalpyinv(sub, heat, state, i)
        # set ∂H∂T (a.k.a ∂H∂T)
        state.∂H∂T[i] = state.T[i] ≈ 0.0 ? 1e8 : state.C[i]
    end
    return nothing
end
