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
        a=Δk[1];
        k*(T₂-T₁)/δ/a
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
        a=Δk[end];
        -k*(T₂-T₁)/δ/a
    end
    return nothing
end

# Default implementation of `variables` for freeze curve
variables(::SubSurface, ::Heat, ::FreezeCurve) = ()
# Fallback (error) implementation for freeze curve
(fc::FreezeCurve)(layer::SubSurface, heat::Heat, state) =
    error("freeze curve $(typeof(fc)) not implemented for $(typeof(heat)) on layer $(typeof(layer))")
freezecurve(heat::Heat) = heat.freezecurve
enthalpy(T, C, L, θ) = T*C + L*θ
heatcapacity!(layer::SubSurface, heat::Heat, state) = error("heatcapacity not defined for $(typeof(heat)) on $(typeof(layer))")
thermalconductivity!(layer::SubSurface, heat::Heat, state) = error("thermalconductivity not defined for $(typeof(heat)) on $(typeof(layer))")
"""
Variable definitions for heat conduction on any subsurface layer. Joins variables defined on
`layer` and `heat` individually as well as variables defined by the freeze curve.
"""
variables(layer::SubSurface, heat::Heat) = (
    variables(layer)..., # layer variables
    variables(heat)...,  # heat variables
    variables(layer, heat, freezecurve(heat))..., # freeze curve variables
)
"""
Variable definitions for heat conduction (enthalpy).
"""
variables(heat::Heat{<:FreezeCurve,Enthalpy}) = (
    Prognostic(:H, Float"J/m^3", OnGrid(Cells)),
    Diagnostic(:T, Float"°C", OnGrid(Cells)),
    Diagnostic(:C, Float"J//K/m^3", OnGrid(Cells)),
    Diagnostic(:Ceff, Float"J/K/m^3", OnGrid(Cells)),
    Diagnostic(:k, Float"W/m/K", OnGrid(Edges)),
    Diagnostic(:kc, Float"W//m/K", OnGrid(Cells)),
    Diagnostic(:θl, Float"1", OnGrid(Cells)),
)
"""
Variable definitions for heat conduction (temperature).
"""
variables(heat::Heat{<:FreezeCurve,Temperature}) = (
    Prognostic(:T, Float"°C", OnGrid(Cells)),
    Diagnostic(:H, Float"J/m^3", OnGrid(Cells)),
    Diagnostic(:dH, Float"J/s/m^3", OnGrid(Cells)),
    Diagnostic(:C, Float"J/K/m^3", OnGrid(Cells)),
    Diagnostic(:Ceff, Float"J/K/m^3", OnGrid(Cells)),
    Diagnostic(:k, Float"W/m/K", OnGrid(Edges)),
    Diagnostic(:kc, Float"W/m/K", OnGrid(Cells)),
    Diagnostic(:θl, Float"1", OnGrid(Cells)),
)
""" Diagonstic step for heat conduction (all state configurations) on any subsurface layer. """
function diagnosticstep!(layer::SubSurface, heat::Heat, state)
    # Reset energy flux to zero; this is redundant when H is the prognostic variable
    # but necessary when it is not.
    @. state.dH = zero(eltype(state.dH))
    # Evaluate the freeze curve (updates T, C, and θl)
    fc! = freezecurve(heat);
    fc!(layer, heat, state)
    # Update thermal conductivity
    thermalconductivity!(layer, heat, state)
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
    heatconduction!(state.dH,state.T,ΔT,state.k,Δk)
end
""" Prognostic step for heat conduction (temperature) on subsurface layer. """
function prognosticstep!(::SubSurface, ::Heat{<:FreezeCurve,Temperature}, state)
    Δk = Δ(state.grids.k) # cell sizes
    ΔT = Δ(state.grids.T)
    # Diffusion on non-boundary cells
    heatconduction!(state.dH,state.T,ΔT,state.k,Δk)
    # Compute temperature flux by dividing by C_eff;
    # C_eff should be computed by the freeze curve.
    @inbounds @. state.dT = state.dH / state.Ceff
    return nothing
end
boundaryflux(::Neumann, bc::BoundaryProcess, top::Top, heat::Heat, sub::SubSurface, stop, ssub) = boundaryvalue(bc,top,heat,sub,stop,ssub)
boundaryflux(::Neumann, bc::BoundaryProcess, bot::Bottom, heat::Heat, sub::SubSurface, sbot, ssub) = boundaryvalue(bc,bot,heat,sub,sbot,ssub)
function boundaryflux(::Dirichlet, bc::BoundaryProcess, top::Top, heat::Heat, sub::SubSurface, stop, ssub)
    Δk = Δ(ssub.grids.k)
    @inbounds let Tupper=boundaryvalue(bc,top,heat,sub,stop,ssub),
        Tsub=ssub.T[1],
        k=ssub.k[1],
        d=Δk[1],
        δ₀=(d/2); # distance to boundary
        -k*(Tsub-Tupper)/δ₀/d
    end
end
function boundaryflux(::Dirichlet, bc::BoundaryProcess, bot::Bottom, heat::Heat, sub::SubSurface, sbot, ssub)
    Δk = Δ(ssub.grids.k)
    @inbounds let Tlower=boundaryvalue(bc,bot,heat,sub,sbot,ssub),
        Tsub=ssub.T[end],
        k=ssub.k[end],
        d=Δk[1],
        δ₀=(d/2); # distance to boundary
        -k*(Tsub-Tlower)/δ₀/d
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
    @inbounds ssub.dH[1] += boundaryflux(bc, top, heat, sub, stop, ssub)
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
    @inbounds ssub.dH[end] += boundaryflux(bc, bot, heat, sub, sbot, ssub)
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
"""
Implementation of "free water" freeze curve for any subsurface layer. Assumes that
'state' contains at least temperature (T), enthalpy (H), heat capacity (C),
total water content (θw), and liquid water content (θl).
"""
@inline function (fc::FreeWater)(layer::SubSurface, heat::Heat{FreeWater,Enthalpy}, state)
    @inline function enthalpyinv(H, C, L, θtot)
        let θtot = max(1.0e-8, θtot),
            Lθ = L*θtot,
            I_t = H > Lθ,
            I_f = H <= 0.0;
            (I_t*(H-Lθ) + I_f*H)/C
        end
    end
    @inline function freezethaw(H, L, θtot)
        let θtot = max(1.0e-8, θtot),
            Lθ = L*θtot,
            I_t = H > Lθ,
            I_c = (H > 0.0) && (H <= Lθ);
            I_c*(H/Lθ) + I_t
        end
    end
    @inbounds @. state.θl = freezethaw(state.H, heat.L, state.θw)*state.θw
    heatcapacity!(layer, heat, state) # update heat capacity, C
    @inbounds @. state.T = enthalpyinv(state.H, state.C, heat.L, state.θw)
    return nothing
end
