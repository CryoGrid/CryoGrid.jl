function heatflux(T₁, T₂, k₁, k₂, Δ₁, Δ₂, z₁, z₂)
    # thermal conductivity between cells
    k = Numerics.harmonicmean(k₁, k₂, Δ₁, Δ₂)
    # calculate heat flux between cells
    jH = @inbounds let δ = z₂ - z₁;
        Numerics.flux(T₁, T₂, δ, k)
    end
    return (; jH, k)
end

# Free water freeze curve
@inline function enthalpyinv(sub::SubSurface, heat::HeatBalance{FreeWater,<:Enthalpy}, state, i)
    θwi = Hydrology.watercontent(sub, state, i)
    H = state.H[i]
    L = heat.prop.L
    θw, I_t, I_f, I_c, Lθ = FreezeCurves.freewater(H, θwi, L)
    θfracs = volumetricfractions(sub, state, i)
    C = heatcapacity(sub, heat, θfracs...)
    T = (I_t*(H-Lθ) + I_f*H)/C
    return T, θw, C
end

"""
    freezethaw!(sub::SubSurface, heat::HeatBalance{FreeWater,<:Enthalpy}, state)

Implementation of "free water" freezing characteristic for any subsurface layer.
Assumes that `state` contains at least temperature (T), enthalpy (H), heat capacity (C),
total water content (θwi), and liquid water content (θw).
"""
freezethaw!(sub, state) = freezethaw!(sub, processes(sub), state)
@inline function freezethaw!(sub::SubSurface, heat::HeatBalance{FreeWater,<:Enthalpy}, state)
    @inbounds for i in 1:length(state.H)
        # update T, θw, C
        state.T[i], state.θw[i], state.C[i] = enthalpyinv(sub, heat, state, i)
        # set ∂H∂T (a.k.a ∂H∂T)
        state.∂H∂T[i] = state.T[i] ≈ 0.0 ? 1e8 : state.C[i]
    end
    return nothing
end

"""
Common variable definitions for all heat implementations.
"""
heatvariables(::HeatBalance) = (
    Diagnostic(:jH, OnGrid(Edges), u"W/m^2"),    
    Diagnostic(:∂H∂T, OnGrid(Cells), u"J/K/m^3", domain=0..Inf),
    Diagnostic(:∂θw∂T, OnGrid(Cells), domain=0..Inf),
    Diagnostic(:C, OnGrid(Cells), u"J/K/m^3"),
    Diagnostic(:k, OnGrid(Edges), u"W/m/K"),
    Diagnostic(:kc, OnGrid(Cells), u"W/m/K"),
    Diagnostic(:θw, OnGrid(Cells), domain=0..1),
)

"""
Variable definitions for heat conduction (enthalpy) on any SubSurface layer.
"""
CryoGrid.variables(heat::HeatBalance{<:FreezeCurve,<:Enthalpy}) = (
    Prognostic(:H, OnGrid(Cells), u"J/m^3"),
    Diagnostic(:T, OnGrid(Cells), u"°C"),
    heatvariables(heat)...,
)
"""
Variable definitions for heat conduction (temperature).
"""
CryoGrid.variables(heat::HeatBalance{<:FreezeCurve,<:Temperature}) = (
    Prognostic(:T, OnGrid(Cells), u"°C"),
    Diagnostic(:H, OnGrid(Cells), u"J/m^3"),
    Diagnostic(:∂H∂t, OnGrid(Cells), u"W/m^3"),
    heatvariables(heat)...,
)

"""
Diagnostic state update for heat conduction (all state configurations) on any subsurface layer.
"""
function CryoGrid.updatestate!(sub::SubSurface, heat::HeatBalance, state)
    resetfluxes!(sub, heat, state)
    # Evaluate freeze/thaw processes
    freezethaw!(sub, state)
    # Update thermal conductivity
    thermalconductivity!(sub, heat, state)
    return nothing
end

function CryoGrid.interact!(sub1::SubSurface, ::HeatBalance, sub2::SubSurface, ::HeatBalance, s1, s2)
    Δk₁ = CryoGrid.thickness(sub1, s1, last)
    Δk₂ = CryoGrid.thickness(sub2, s2, first)
    z₁ = CryoGrid.midpoint(sub1, s1, last)
    z₂ = CryoGrid.midpoint(sub2, s2, first)
    jH, k = heatflux(s1.T[end], s2.T[1], s1.kc[end], s2.kc[1], Δk₁, Δk₂, z₁, z₂)
    # set edge conductivity
    s2.k[1] = s1.k[end] = k
    # set inner boundary flux
    s2.jH[1] = s1.jH[end] += jH
    return nothing
end


function CryoGrid.computefluxes!(::SubSurface, ::HeatBalance{<:FreezeCurve,<:Enthalpy}, state)
    Δk = Δ(state.grid) # cell sizes
    ΔT = Δ(cells(state.grid)) # midpoint distances
    # compute internal fluxes and non-linear diffusion assuming boundary fluxes have been set
    Numerics.nonlineardiffusion!(state.∂H∂t, state.jH, state.T, ΔT, state.k, Δk)
    return nothing
end
function CryoGrid.computefluxes!(sub::SubSurface, ::HeatBalance{<:FreezeCurve,<:Temperature}, state)
    Δk = Δ(state.grid) # cell sizes
    ΔT = Δ(cells(state.grid)) # midpoint distances
    # compute internal fluxes and non-linear diffusion assuming boundary fluxes have been set
    Numerics.nonlineardiffusion!(state.∂H∂t, state.jH, state.T, ΔT, state.k, Δk)
    # Compute temperature flux by dividing by ∂H∂T;
    # ∂H∂T should be computed by the freeze curve.
    @inbounds @. state.∂T∂t = state.∂H∂t / state.∂H∂T
    return nothing
end

# CFL not defined for free-water freeze curve
CryoGrid.timestep(::SubSurface, heat::HeatBalance{FreeWater,<:Enthalpy,<:CryoGrid.CFL}, state) = Inf
"""
    timestep(::SubSurface, ::HeatBalance{Tfc,THeatOp,CFL}, state) where {Tfc,THeatOp}

Implementation of `timestep` for `HeatBalance` using the Courant-Fredrichs-Lewy condition
defined as: Δt_max = u*Δx^2, where`u` is the "characteristic velocity" which here
is taken to be the diffusivity: `∂H∂T / kc`.
"""
function CryoGrid.timestep(::SubSurface, heat::HeatBalance{Tfc,THeatOp,<:CryoGrid.CFL}, state) where {Tfc,THeatOp}
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
function CryoGrid.timestep(::SubSurface, heat::HeatBalance{Tfc,THeatOp,<:CryoGrid.MaxDelta}, state) where {Tfc,THeatOp}
    Δx = Δ(state.grid)
    dtmax = Inf
    @inbounds for i in eachindex(Δx)
        dtmax = min(dtmax, heat.dtlim(state.∂H∂t[i], state.H[i], state.t))
    end
    dtmax = isfinite(dtmax) && dtmax > 0 ? dtmax : Inf
    return dtmax
end

function CryoGrid.resetfluxes!(sub::SubSurface, heat::HeatBalance, state)
    # Reset energy fluxes to zero; this is redundant when H is the prognostic variable
    # but necessary when it is not.
    @. state.∂H∂t = zero(eltype(state.∂H∂t))
    @. state.jH = zero(eltype(state.jH))
end
