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
    totalwater(sub::SubSurface, state)
    totalwater(sub::SubSurface, heat::Heat, state)
    totalwater(sub::SubSurface, heat::Heat, state, i)

Retrieves the total water content for the given layer at grid cell `i`, if specified.
Defaults to using the scalar total water content defined on layer `sub`.
"""
@inline totalwater(sub::SubSurface, state) = error("totalwater not defined for $(typeof(sub))")
@inline totalwater(sub::SubSurface, heat::Heat, state) = totalwater(sub, state)
@inline totalwater(sub::SubSurface, heat::Heat, state, i) = Utils.getscalar(totalwater(sub, heat, state), i)
"""
    liquidwater(sub::SubSurface, heat::Heat, state)
    liquidwater(sub::SubSurface, heat::Heat, state, i)

Retrieves the liquid water content for the given layer at grid cell `i`, if specified.
Defaults to using the scalar total water content defined on layer `sub`.
"""
@inline liquidwater(sub::SubSurface, heat::Heat, state) = state.θl
@inline liquidwater(sub::SubSurface, heat::Heat, state, i) = Utils.getscalar(liquidwater(sub, heat, state), i)
"""
    heatcapacity(capacities::NTuple{N}, fracs::NTuple{N}) where N

Computes the heat capacity as a weighted average over constituent `capacities` with volumetric fractions `fracs`.
"""
heatcapacity(capacities::NTuple{N}, fracs::NTuple{N}) where N = sum(map(*, capacities, fracs))
"""
    heatcapacity!(sub::SubSurface, heat::Heat, state)

Computes the heat capacity for the given layer from the current state and stores the result in-place in the state variable `C`.
"""
@inline function heatcapacity!(sub::SubSurface, heat::Heat, state)
    @inbounds for i in 1:length(state.T)
        state.C[i] = heatcapacity(sub, heat, state, i)
    end
end
"""
    thermalconductivity(capacities::NTuple{N}, fracs::NTuple{N}) where N

Computes the thermal conductivity as a squared weighted sum over constituent `conductivities` with volumetric fractions `fracs`.
"""
thermalconductivity(conductivities::NTuple{N}, fracs::NTuple{N}) where N = sum(map(*, map(sqrt, conductivities), fracs))^2
"""
    thermalconductivity!(sub::SubSurface, heat::Heat, state)

Computes the thermal conductivity for the given layer from the current state and stores the result in-place in the state variable `C`.
"""
@inline function thermalconductivity!(sub::SubSurface, heat::Heat, state)
    @inbounds for i in 1:length(state.T)
        state.kc[i] = thermalconductivity(sub, heat, state, i)
    end
end
"""
    heatconduction!(∂H,T,ΔT,k,Δk)

1-D heat conduction/diffusion given T, k, and their deltas. Resulting enthalpy gradient is stored in ∂H.
Note that this function does not perform bounds checking. It is up to the user to ensure that all variables are
arrays of the correct length.
"""
function heatconduction!(∂H,T,ΔT,k,Δk)
    # upper boundary
    @inbounds ∂H[1] += let T₂=T[2],
        T₁=T[1],
        k=k[2],
        δ=ΔT[1],
        d=Δk[1];
        k*(T₂-T₁)/δ/d
    end
    # diffusion on non-boundary cells
    @inbounds let T = T,
        k = (@view k[2:end-1]),
        Δk = (@view Δk[2:end-1]),
        ∂H = (@view ∂H[2:end-1]);
        nonlineardiffusion!(∂H, T, ΔT, k, Δk)
    end
    # lower boundary
    @inbounds ∂H[end] += let T₂=T[end],
        T₁=T[end-1],
        k=k[end-1],
        δ=ΔT[end],
        d=Δk[end];
        -k*(T₂-T₁)/δ/d
    end
    return nothing
end
"""
Variable definitions for heat conduction on any subsurface layer. Joins variables defined on
`layer` and `heat` individually as well as variables defined by the freeze curve.
"""
variables(sub::SubSurface, heat::Heat) = (
    variables(sub)..., # layer variables
    variables(heat)...,  # heat variables
    variables(sub, heat, freezecurve(heat))..., # freeze curve variables
)
"""
Variable definitions for heat conduction (enthalpy).
"""
variables(heat::Heat{<:FreezeCurve,Enthalpy}) = (
    Prognostic(:H, Float"J/m^3", OnGrid(Cells)),
    Diagnostic(:T, Float"°C", OnGrid(Cells)),
    basevariables(heat)...,
)
"""
Variable definitions for heat conduction (temperature).
"""
variables(heat::Heat{<:FreezeCurve,Temperature}) = (
    Prognostic(:T, Float"°C", OnGrid(Cells)),
    Diagnostic(:H, Float"J/m^3", OnGrid(Cells)),
    Diagnostic(:dH, Float"W/m^3", OnGrid(Cells)),
    basevariables(heat)...,
)
"""
Common variable definitions for all heat implementations.
"""
basevariables(::Heat) = (
    Diagnostic(:dHdT, Float"J/K/m^3", OnGrid(Cells)),
    Diagnostic(:C, Float"J/K/m^3", OnGrid(Cells)),
    Diagnostic(:k, Float"W/m/K", OnGrid(Edges)),
    Diagnostic(:kc, Float"W/m/K", OnGrid(Cells)),
    Diagnostic(:θl, Float"1", OnGrid(Cells)),
)
""" Diagonstic step for heat conduction (all state configurations) on any subsurface layer. """
function diagnosticstep!(sub::SubSurface, heat::Heat, state)
    # Reset energy flux to zero; this is redundant when H is the prognostic variable
    # but necessary when it is not.
    @. state.dH = zero(eltype(state.dH))
    # Evaluate the freeze curve (updates T, C, and θl)
    fc! = freezecurve(heat);
    fc!(sub, heat, state)
    # Update thermal conductivity
    thermalconductivity!(sub, heat, state)
    # Harmonic mean of inner conductivities
    @inbounds let k = (@view state.k[2:end-1]),
        Δk = Δ(state.grids.k);
        harmonicmean!(k, state.kc, Δk)
    end
    return nothing # ensure no allocation
end
""" Prognostic step for heat conduction (enthalpy) on subsurface layer. """
function prognosticstep!(::SubSurface, ::Heat{<:FreezeCurve,Enthalpy}, state)
    Δk = Δ(state.grids.k) # cell sizes
    ΔT = Δ(state.grids.T)
    # Diffusion on non-boundary cells
    heatconduction!(state.dH, state.T, ΔT, state.k, Δk)
end
""" Prognostic step for heat conduction (temperature) on subsurface layer. """
function prognosticstep!(::SubSurface, ::Heat{<:FreezeCurve,Temperature}, state)
    Δk = Δ(state.grids.k) # cell sizes
    ΔT = Δ(state.grids.T)
    # Diffusion on non-boundary cells
    heatconduction!(state.dH, state.T, ΔT, state.k, Δk)
    # Compute temperature flux by dividing by dHdT;
    # dHdT should be computed by the freeze curve.
    @inbounds @. state.dT = state.dH / state.dHdT
    return nothing
end
@inline boundaryflux(::Neumann, bc::BoundaryProcess, top::Top, heat::Heat, sub::SubSurface, stop, ssub) = boundaryvalue(bc,top,heat,sub,stop,ssub)
@inline boundaryflux(::Neumann, bc::BoundaryProcess, bot::Bottom, heat::Heat, sub::SubSurface, sbot, ssub) = boundaryvalue(bc,bot,heat,sub,sbot,ssub)
@inline function boundaryflux(::Dirichlet, bc::BoundaryProcess, top::Top, heat::Heat, sub::SubSurface, stop, ssub)
    Δk = Δ(ssub.grids.k)
    @inbounds let Tupper=boundaryvalue(bc,top,heat,sub,stop,ssub),
        Tsub=ssub.T[1],
        k=ssub.k[1],
        d=Δk[1],
        δ₀=d/2; # distance to boundary
        -k*(Tsub-Tupper)/δ₀
    end
end
@inline function boundaryflux(::Dirichlet, bc::BoundaryProcess, bot::Bottom, heat::Heat, sub::SubSurface, sbot, ssub)
    Δk = Δ(ssub.grids.k)
    @inbounds let Tlower=boundaryvalue(bc,bot,heat,sub,sbot,ssub),
        Tsub=ssub.T[end],
        k=ssub.k[end],
        d=Δk[end],
        δ₀=d/2; # distance to boundary
        -k*(Tsub-Tlower)/δ₀
    end
end
"""
Generic top interaction. Computes flux dH at top cell.
"""
function interact!(top::Top, bc::BoundaryProcess, sub::SubSurface, heat::Heat, stop, ssub)
    # thermal conductivity at boundary
    # assumes (1) k has already been computed, (2) surface conductivity = cell conductivity
    @inbounds ssub.k[1] = ssub.kc[1]
    # boundary flux
    @log dH_upper = boundaryflux(bc, top, heat, sub, stop, ssub)
    Δk = Δ(ssub.grids.k)
    @inbounds ssub.dH[1] += dH_upper / Δk[1]
    return nothing # ensure no allocation
end
"""
Generic bottom interaction. Computes flux dH at bottom cell.
"""
function interact!(sub::SubSurface, heat::Heat, bot::Bottom, bc::BoundaryProcess, ssub, sbot)
    # thermal conductivity at boundary
    # assumes (1) k has already been computed, (2) bottom conductivity = cell conductivity
    @inbounds ssub.k[end] = ssub.kc[end]
    # boundary flux
    @log dH_lower = boundaryflux(bc, bot, heat, sub, sbot, ssub)
    Δk = Δ(ssub.grids.k)
    @inbounds ssub.dH[end] += dH_lower / Δk[end]
    return nothing # ensure no allocation
end
"""
Generic subsurface interaction. Computes flux dH at boundary between subsurface layers.
"""
function interact!(::SubSurface, ::Heat, ::SubSurface, ::Heat, s1, s2)
    # thermal conductivity between cells
    @inbounds let k₁ = s1.kc[end],
        k₂ = s2.kc[1],
        Δ₁ = Δ(s1.grids.k)[end],
        Δ₂ = Δ(s2.grids.k)[1];
        k = harmonicmean(k₁, k₂, Δ₁, Δ₂);
        s1.k[end] = s2.k[1] = k
    end
    # calculate heat flux between cells
    Qᵢ = @inbounds let k₁ = s1.k[end],
        k₂ = s2.k[1],
        k = k₁ = k₂,
        δ = s2.grids.T[1] - s1.grids.T[end];
        k*(s2.T[1] - s1.T[end]) / δ
    end
    # add fluxes scaled by grid cell size
    @inbounds s1.dH[end] += Qᵢ / Δ(s1.grids.k)[end]
    @inbounds s2.dH[1] += -Qᵢ / Δ(s2.grids.k)[1]
    return nothing # ensure no allocation
end
# Free water freeze curve
@inline function enthalpyinv(sub::SubSurface, heat::Heat{FreeWater,Enthalpy}, state, i)
    let θtot = totalwater(sub, heat, state, i),
        H = state.H[i],
        C = state.C[i],
        L = heat.L,
        Lθ = L*θtot,
        I_t = H > Lθ,
        I_f = H <= 0.0;
        (I_t*(H-Lθ) + I_f*H)/C
    end
end
@inline function liquidwater(sub::SubSurface, heat::Heat{FreeWater,Enthalpy}, state, i)
    let θtot = max(1e-8, totalwater(sub, heat, state, i)),
        H = state.H[i],
        L = heat.L,
        Lθ = L*θtot,
        I_t = H > Lθ,
        I_c = (H > 0.0) && (H <= Lθ);
        (I_c*(H/Lθ) + I_t)θtot
    end
end
"""
    (fc::FreeWater)(sub::SubSurface, heat::Heat{FreeWater,Enthalpy}, state)

Implementation of "free water" freeze curve for any subsurface layer. Assumes that
'state' contains at least temperature (T), enthalpy (H), heat capacity (C),
total water content (θw), and liquid water content (θl).
"""
@inline function (fc::FreeWater)(sub::SubSurface, heat::Heat{FreeWater,Enthalpy}, state)
    @inbounds for i in 1:length(state.H)
        # liquid water content = (total water content) * (liquid fraction)
        state.θl[i] = liquidwater(sub, heat, state, i)
        # update heat capacity
        state.C[i] = heatcapacity(sub, heat, state, i)
        # enthalpy inverse function
        state.T[i] = enthalpyinv(sub, heat, state, i)
        # set dHdT (a.k.a dHdT)
        state.dHdT[i] = state.T[i] ≈ 0.0 ? 1e8 : 1/state.C[i]
    end
    return nothing
end
