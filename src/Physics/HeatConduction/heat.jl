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
    C_eff(T, C, L, dθdT, cw, ci) = C + dθdT*(L + T*(cw - ci))

Computes the apparent or "effective" heat capacity `dHdT` as a function of temperature, volumetric heat capacity,
latent heat of fusion, derivative of the freeze curve `dθdT`, and the constituent heat capacities of water and ice.
"""
@inline C_eff(T, C, L, dθdT, cw, ci) = C + dθdT*(L + T*(cw - ci))
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

# Generic heat conduction implementation

"""
    volumetricfractions(::SubSurface, heat::Heat, state, i)

Get constituent volumetric fractions for generic `SubSurface` layer. Default implementation assumes
the only constitutents are liquid water, ice, and air: `(θw,θi,θa)`.
"""
@inline function Physics.volumetricfractions(sub::SubSurface, heat::Heat, state, i)
    let θwi = waterice(sub, heat, state, i),
        θw = liquidwater(sub, heat, state, i),
        θa = 1.0 - θwi,
        θi = θwi - θw;
        return (θw, θi, θa)
    end
end
"""
    freezethaw!(sub::SubSurface, heat::Heat, state)

Calculates freezing and thawing effects, including evaluation of the freeze curve.
In general, this function should compute at least the liquid/frozen water contents
and the corresponding heat capacity. Other variables such as temperature or enthalpy
should also be computed depending on the thermal scheme being implemented.
"""
freezethaw!(sub::SubSurface, heat::Heat, state) = error("missing implementation of freezethaw!")
"""
    heatcapacity(sub::SubSurface, heat::Heat, θfracs...)

Computes the heat capacity as a weighted average over constituent capacities with volumetric fractions `θfracs`.
"""
@inline function heatcapacity(sub::SubSurface, heat::Heat, θfracs...)
    cs = heatcapacities(sub, heat)
    return sum(map(*, cs, θfracs))
end
"""
    heatcapacity!(sub::SubSurface, heat::Heat, state)

Computes the heat capacity for the given layer from the current state and stores the result in-place in the state variable `C`.
"""
@inline function heatcapacity!(sub::SubSurface, heat::Heat, state)
    @inbounds for i in 1:length(state.T)
        θfracs = volumetricfractions(sub, heat, state, i)
        state.C[i] = heatcapacity(sub, heat, θfracs...)
    end
end
"""
    thermalconductivity(sub::SubSurface, heat::Heat, θfracs...)

Computes the thermal conductivity as a squared weighted sum over constituent conductivities with volumetric fractions `θfracs`.
"""
@inline function thermalconductivity(sub::SubSurface, heat::Heat, θfracs...)
    ks = thermalconductivities(sub, heat)
    return sum(map(*, map(sqrt, ks), θfracs))^2
end
"""
    thermalconductivity!(sub::SubSurface, heat::Heat, state)

Computes the thermal conductivity for the given layer from the current state and stores the result in-place in the state variable `C`.
"""
@inline function thermalconductivity!(sub::SubSurface, heat::Heat, state)
    @inbounds for i in 1:length(state.T)
        θfracs = volumetricfractions(sub, heat, state, i)
        state.kc[i] = thermalconductivity(sub, heat, θfracs...)
    end
end
"""
Variable definitions for heat conduction on any subsurface layer. Joins variables defined on
`layer` and `heat` individually as well as variables defined by the freeze curve.
"""
CryoGrid.variables(sub::SubSurface, heat::Heat) = (
    variables(sub)..., # layer variables
    variables(heat)...,  # heat variables
    variables(sub, heat, freezecurve(heat))..., # freeze curve variables
)
"""
Variable definitions for heat conduction (enthalpy).
"""
CryoGrid.variables(heat::Heat{<:FreezeCurve,Enthalpy}) = (
    Prognostic(:H, OnGrid(Cells), u"J/m^3"),
    Diagnostic(:T, OnGrid(Cells), u"°C"),
    heatvariables(heat)...,
)
"""
Variable definitions for heat conduction (temperature).
"""
CryoGrid.variables(heat::Heat{<:FreezeCurve,Temperature}) = (
    Prognostic(:T, OnGrid(Cells), u"°C"),
    Diagnostic(:H, OnGrid(Cells), u"J/m^3"),
    Diagnostic(:dH, OnGrid(Cells), u"W/m^3"),
    Diagnostic(:dθdT, OnGrid(Cells)),
    heatvariables(heat)...,
)
"""
Common variable definitions for all heat implementations.
"""
heatvariables(::Heat) = (
    Diagnostic(:dH_upper, Scalar, u"J/K/m^2"),
    Diagnostic(:dH_lower, Scalar, u"J/K/m^2"),
    Diagnostic(:dHdT, OnGrid(Cells), u"J/K/m^3"),
    Diagnostic(:C, OnGrid(Cells), u"J/K/m^3"),
    Diagnostic(:k, OnGrid(Edges), u"W/m/K"),
    Diagnostic(:kc, OnGrid(Cells), u"W/m/K"),
    Diagnostic(:θw, OnGrid(Cells)),
)
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
Diagonstic step for heat conduction (all state configurations) on any subsurface layer.
"""
function CryoGrid.diagnosticstep!(sub::SubSurface, heat::Heat, state)
    # Reset energy flux to zero; this is redundant when H is the prognostic variable
    # but necessary when it is not.
    @. state.dH = zero(eltype(state.dH))
    @. state.dH_upper = zero(eltype(state.dH_upper))
    @. state.dH_lower = zero(eltype(state.dH_lower))
    # Evaluate freeze/thaw processes
    freezethaw!(sub, heat, state)
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
"""
Generic top interaction. Computes flux dH at top cell.
"""
function CryoGrid.interact!(top::Top, bc::HeatBC, sub::SubSurface, heat::Heat, stop, ssub)
    Δk = CryoGrid.thickness(sub, ssub, first) # using `thickness` allows for generic layer implementations
    # boundary flux
    @setscalar ssub.dH_upper = boundaryflux(bc, top, heat, sub, stop, ssub)
    @inbounds ssub.dH[1] += getscalar(ssub.dH_upper) / Δk
    return nothing # ensure no allocation
end
"""
Generic bottom interaction. Computes flux dH at bottom cell.
"""
function CryoGrid.interact!(sub::SubSurface, heat::Heat, bot::Bottom, bc::HeatBC, ssub, sbot)
    Δk = CryoGrid.thickness(sub, ssub, last) # using `thickness` allows for generic layer implementations
    # boundary flux
    @setscalar ssub.dH_lower = boundaryflux(bc, bot, heat, sub, sbot, ssub)
    @inbounds ssub.dH[end] += getscalar(ssub.dH_lower) / Δk
    return nothing # ensure no allocation
end
"""
Generic subsurface interaction. Computes flux dH at boundary between subsurface layers.
"""
function CryoGrid.interact!(sub1::SubSurface, ::Heat, sub2::SubSurface, ::Heat, s1, s2)
    Δk₁ = CryoGrid.thickness(sub1, s1, first)
    Δk₂ = CryoGrid.thickness(sub2, s2, last)
    # thermal conductivity between cells
    k = s1.k[end] = s2.k[1] =
        @inbounds let k₁ = s1.kc[end],
            k₂ = s2.kc[1],
            Δ₁ = Δk₁[end],
            Δ₂ = Δk₂[1];
            harmonicmean(k₁, k₂, Δ₁, Δ₂)
        end
    # calculate heat flux between cells
    Qᵢ = @inbounds let z₁ = CryoGrid.midpoint(sub1, s1, last),
        z₂ = CryoGrid.midpoint(sub2, s2, first),
        δ = z₂ - z₁;
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
"""
Prognostic step for heat conduction (enthalpy) on subsurface layer.
"""
function CryoGrid.prognosticstep!(::SubSurface, ::Heat{<:FreezeCurve,Enthalpy}, state)
    Δk = Δ(state.grids.k) # cell sizes
    ΔT = Δ(state.grids.T)
    # Diffusion on non-boundary cells
    heatconduction!(state.dH, state.T, ΔT, state.k, Δk)
end
"""
Prognostic step for heat conduction (temperature) on subsurface layer.
"""
function CryoGrid.prognosticstep!(sub::SubSurface, ::Heat{<:FreezeCurve,Temperature}, state)
    Δk = Δ(state.grids.k) # cell sizes
    ΔT = Δ(state.grids.T)
    # Diffusion on non-boundary cells
    heatconduction!(state.dH, state.T, ΔT, state.k, Δk)
    # Compute temperature flux by dividing by dHdT;
    # dHdT should be computed by the freeze curve.
    @inbounds @. state.dT = state.dH / state.dHdT
    return nothing
end
# Boundary fluxes
@inline CryoGrid.boundaryflux(::Neumann, bc::HeatBC, top::Top, heat::Heat, sub::SubSurface, stop, ssub) = boundaryvalue(bc,top,heat,sub,stop,ssub)
@inline CryoGrid.boundaryflux(::Neumann, bc::HeatBC, bot::Bottom, heat::Heat, sub::SubSurface, sbot, ssub) = boundaryvalue(bc,bot,heat,sub,sbot,ssub)
@inline function CryoGrid.boundaryflux(::Dirichlet, bc::HeatBC, top::Top, heat::Heat, sub::SubSurface, stop, ssub)
    Δk = CryoGrid.thickness(sub, ssub, first) # using `thickness` allows for generic layer implementations
    @inbounds let Tupper=boundaryvalue(bc,top,heat,sub,stop,ssub),
        Tsub=ssub.T[1],
        k=ssub.k[1],
        δ=Δk/2; # distance to boundary
        -k*(Tsub-Tupper)/δ
    end
end
@inline function CryoGrid.boundaryflux(::Dirichlet, bc::HeatBC, bot::Bottom, heat::Heat, sub::SubSurface, sbot, ssub)
    Δk = CryoGrid.thickness(sub, ssub, last) # using `thickness` allows for generic layer implementations
    @inbounds let Tlower=boundaryvalue(bc,bot,heat,sub,sbot,ssub),
        Tsub=ssub.T[end],
        k=ssub.k[end],
        δ=Δk/2; # distance to boundary
        -k*(Tsub-Tlower)/δ
    end
end
"""
    timestep(::SubSurface, ::Heat{Tfc,TPara,CFL}, state) where {TPara}

Implementation of `timestep` for `Heat` using the Courant-Fredrichs-Lewy condition
defined as: Δt_max = u*Δx^2, where`u` is the "characteristic velocity" which here
is taken to be the diffusivity: `dHdT / kc`.
"""
@inline function CryoGrid.timestep(::SubSurface, heat::Heat{Tfc,TPara,CFL}, state) where {Tfc,TPara}
    Δx = Δ(state.grid)
    dtmax = Inf
    @inbounds for i in 1:length(Δx)
        dtmax = let u = state.dHdT[i] / state.kc[i],
            Δx = Δx[i];
            Δt = Δx^2*u
            min(dtmax, Δt)
        end
    end
    dtmax = isfinite(dtmax) && dtmax > 0 ? dtmax : heat.dtlim.fallback_dt
    return dtmax
end
# Free water freeze curve
@inline function enthalpyinv(sub::SubSurface, heat::Heat{FreeWater,Enthalpy}, state, i)
    f_hc = partial(heatcapacity, liquidwater, sub, heat, state, i)
    return enthalpyinv(heat.freezecurve, f_hc, state.H[i], heat.prop.L, f_hc.θwi)
end
@inline function enthalpyinv(::FreeWater, f_hc::Function, H, L, θtot)
    let θtot = max(1e-8, θtot),
        Lθ = L*θtot,
        I_t = H > Lθ,
        I_f = H <= 0.0,
        I_c = (H > 0.0) && (H <= Lθ);
        # compute liquid water content -> heat capacity -> temperature
        θw = (I_c*(H/Lθ) + I_t)θtot
        C = f_hc(θw)
        T = (I_t*(H-Lθ) + I_f*H)/C
        return T, θw, C
    end
end
"""
    freezethaw!(sub::SubSurface, heat::Heat{FreeWater,Enthalpy}, state)

Implementation of "free water" freezing characteristic for any subsurface layer.
Assumes that `state` contains at least temperature (T), enthalpy (H), heat capacity (C),
total water content (θwi), and liquid water content (θw).
"""
@inline function freezethaw!(sub::SubSurface, heat::Heat{FreeWater,Enthalpy}, state)
    @inbounds for i in 1:length(state.H)
        # update T, θw, C
        state.T[i], state.θw[i], state.C[i] = enthalpyinv(sub, heat, state, i)
        # set dHdT (a.k.a dHdT)
        state.dHdT[i] = state.T[i] ≈ 0.0 ? 1e6 : 1/state.C[i]
    end
    return nothing
end
