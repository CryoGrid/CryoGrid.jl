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
    liquidwater(sub::SubSurface, heat::Heat, state)
    liquidwater(sub::SubSurface, heat::Heat, state, i)

Retrieves the liquid water content for the given layer at grid cell `i`, if specified.
Defaults to using the scalar total water content defined on layer `sub`.
"""
@inline liquidwater(::SubSurface, ::Heat, state) = state.θl
@inline liquidwater(sub::SubSurface, heat::Heat, state, i) = Utils.getscalar(liquidwater(sub, heat, state), i)
"""
    thermalconductivities(::SubSurface, heat::Heat)

Get thermal conductivities for generic `SubSurface` layer.
"""
@inline thermalconductivities(::SubSurface, heat::Heat) = (heat.prop.kw, heat.prop.ki, heat.prop.ka)
"""
    heatcapacities(::SubSurface, heat::Heat)

Get heat capacities for generic `SubSurface` layer.
"""
@inline heatcapacities(::SubSurface, heat::Heat) = (heat.prop.cw, heat.prop.ci, heat.prop.ca)
"""
    volumetricfractions(::SubSurface, heat::Heat)

Get constituent volumetric fractions for generic `SubSurface` layer.
"""
@inline function volumetricfractions(sub::SubSurface, heat::Heat, state, i)
    return let θw = totalwater(sub, state, i),
        θl = liquidwater(sub, heat, state, i),
        θa = 1.0 - θw,
        θi = θw - θl;
        (θl, θi, θa)
    end
end

# Generic heat conduction implementation

"""
    heatcapacity(capacities::NTuple{N}, fracs::NTuple{N}) where N

Computes the heat capacity as a weighted average over constituent `capacities` with volumetric fractions `fracs`.
"""
heatcapacity(capacities::NTuple{N,Any}, fracs::NTuple{N,Any}) where N = sum(map(*, capacities, fracs))
@inline function heatcapacity(sub::SubSurface, heat::Heat, state, i)
    θs = volumetricfractions(sub, heat, state, i)
    cs = heatcapacities(sub, heat)
    return heatcapacity(cs, θs)
end
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
    thermalconductivity(conductivities::NTuple{N}, fracs::NTuple{N}) where N

Computes the thermal conductivity as a squared weighted sum over constituent `conductivities` with volumetric fractions `fracs`.
"""
thermalconductivity(conductivities::NTuple{N,Any}, fracs::NTuple{N,Any}) where N = sum(map(*, map(sqrt, conductivities), fracs))^2
@inline function thermalconductivity(sub::SubSurface, heat::Heat, state, i)
    θs = volumetricfractions(sub, heat, state, i)
    ks = thermalconductivities(sub, heat)
    return thermalconductivity(ks, θs)
end
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
        ϵ=ΔT[1],
        δ=Δk[1];
        k*(T₂-T₁)/ϵ/δ
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
        ϵ=ΔT[end],
        δ=Δk[end];
        -k*(T₂-T₁)/ϵ/δ
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
    Prognostic(:H, OnGrid(Cells), u"J/m^3"),
    Diagnostic(:T, OnGrid(Cells), u"°C"),
    basevariables(heat)...,
)
"""
Variable definitions for heat conduction (temperature).
"""
variables(heat::Heat{<:FreezeCurve,Temperature}) = (
    Prognostic(:T, OnGrid(Cells), u"°C"),
    Diagnostic(:H, OnGrid(Cells), u"J/m^3"),
    Diagnostic(:dH, OnGrid(Cells), u"W/m^3"),
    basevariables(heat)...,
)
"""
Common variable definitions for all heat implementations.
"""
basevariables(::Heat) = (
    Diagnostic(:dH_upper, Scalar, u"J/K/m^2"),
    Diagnostic(:dH_lower, Scalar, u"J/K/m^2"),
    Diagnostic(:dHdT, OnGrid(Cells), u"J/K/m^3"),
    Diagnostic(:C, OnGrid(Cells), u"J/K/m^3"),
    Diagnostic(:k, OnGrid(Edges), u"W/m/K"),
    Diagnostic(:kc, OnGrid(Cells), u"W/m/K"),
    Diagnostic(:θl, OnGrid(Cells)),
)
""" Diagonstic step for heat conduction (all state configurations) on any subsurface layer. """
function diagnosticstep!(sub::SubSurface, heat::Heat, state)
    # Reset energy flux to zero; this is redundant when H is the prognostic variable
    # but necessary when it is not.
    @. state.dH = zero(eltype(state.dH))
    @. state.dH_upper = zero(eltype(state.dH_upper))
    @. state.dH_lower = zero(eltype(state.dH_lower))
    # Evaluate the freeze curve (updates T, C, and θl)
    fc! = freezecurve(heat);
    fc!(sub, heat, state)
    # Update thermal conductivity
    thermalconductivity!(sub, heat, state)
    # thermal conductivity at boundaries
    # assumes boundary conductivities = cell conductivities
    @inbounds state.k[1] = state.kc[1]
    @inbounds state.k[end] = state.kc[end]
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
function prognosticstep!(sub::SubSurface, ::Heat{<:FreezeCurve,Temperature}, state)
    Δk = Δ(state.grids.k) # cell sizes
    ΔT = Δ(state.grids.T)
    # Diffusion on non-boundary cells
    heatconduction!(state.dH, state.T, ΔT, state.k, Δk)
    # Compute temperature flux by dividing by dHdT;
    # dHdT should be computed by the freeze curve.
    @inbounds @. state.dT = state.dH / state.dHdT
    return nothing
end
@inline boundaryflux(::Neumann, bc::HeatBC, top::Top, heat::Heat, sub::SubSurface, stop, ssub) = boundaryvalue(bc,top,heat,sub,stop,ssub)
@inline boundaryflux(::Neumann, bc::HeatBC, bot::Bottom, heat::Heat, sub::SubSurface, sbot, ssub) = boundaryvalue(bc,bot,heat,sub,sbot,ssub)
@inline function boundaryflux(::Dirichlet, bc::HeatBC, top::Top, heat::Heat, sub::SubSurface, stop, ssub)
    Δk = thickness(sub, ssub) # using `thickness` allows for generic layer implementations
    @inbounds let Tupper=boundaryvalue(bc,top,heat,sub,stop,ssub),
        Tsub=ssub.T[1],
        k=ssub.k[1],
        δ=Δk[1]/2; # distance to boundary
        -k*(Tsub-Tupper)/δ
    end
end
@inline function boundaryflux(::Dirichlet, bc::HeatBC, bot::Bottom, heat::Heat, sub::SubSurface, sbot, ssub)
    Δk = thickness(sub, ssub) # using `thickness` allows for generic layer implementations
    @inbounds let Tlower=boundaryvalue(bc,bot,heat,sub,sbot,ssub),
        Tsub=ssub.T[end],
        k=ssub.k[end],
        δ=Δk[end]/2; # distance to boundary
        -k*(Tsub-Tlower)/δ
    end
end
"""
Generic top interaction. Computes flux dH at top cell.
"""
function interact!(top::Top, bc::HeatBC, sub::SubSurface, heat::Heat, stop, ssub)
    Δk = thickness(sub, ssub) # using `thickness` allows for generic layer implementations
    # boundary flux
    @setscalar ssub.dH_upper = boundaryflux(bc, top, heat, sub, stop, ssub)
    @inbounds ssub.dH[1] += getscalar(ssub.dH_upper) / Δk[1]
    return nothing # ensure no allocation
end
"""
Generic bottom interaction. Computes flux dH at bottom cell.
"""
function interact!(sub::SubSurface, heat::Heat, bot::Bottom, bc::HeatBC, ssub, sbot)
    Δk = thickness(sub, ssub) # using `thickness` allows for generic layer implementations
    # boundary flux
    @setscalar ssub.dH_lower = boundaryflux(bc, bot, heat, sub, sbot, ssub)
    @inbounds ssub.dH[end] += getscalar(ssub.dH_lower) / Δk[end]
    return nothing # ensure no allocation
end
"""
Generic subsurface interaction. Computes flux dH at boundary between subsurface layers.
"""
function interact!(sub1::SubSurface, ::Heat, sub2::SubSurface, ::Heat, s1, s2)
    Δk₁ = thickness(sub1, s1)
    Δk₂ = thickness(sub2, s2)
    # thermal conductivity between cells
    @inbounds let k₁ = s1.kc[end],
        k₂ = s2.kc[1],
        Δ₁ = Δk₁[end],
        Δ₂ = Δk₂[1];
        k = harmonicmean(k₁, k₂, Δ₁, Δ₂);
        s1.k[end] = s2.k[1] = k
    end
    # calculate heat flux between cells
    Qᵢ = @inbounds let k₁ = s1.k[end],
        k₂ = s2.k[1],
        k = k₁ = k₂,
        z₁ = midpoints(sub1, s1),
        z₂ = midpoints(sub2, s2),
        δ = z₂[1] - z₁[end];
        k*(s2.T[1] - s1.T[end]) / δ
    end
    # diagnostics
    @setscalar s1.dH_lower = Qᵢ
    @setscalar s2.dH_upper = -Qᵢ
    # add fluxes scaled by grid cell size
    @inbounds s1.dH[end] += Qᵢ / Δk₁[end]
    @inbounds s2.dH[1] += -Qᵢ / Δk₂[1]
    return nothing # ensure no allocation
end
# Free water freeze curve
@inline function enthalpyinv(sub::SubSurface, heat::Heat{FreeWater,Enthalpy}, state, i)
    let θtot = totalwater(sub, state, i),
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
    let θtot = max(1e-8, totalwater(sub, state, i)),
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
