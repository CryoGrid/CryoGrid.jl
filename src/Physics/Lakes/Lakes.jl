module Lakes

import CryoGrid

using CryoGrid
using CryoGrid.Physics.Heat
using CryoGrid.Utils
using CryoGrid.Numerics

using Unitful

export Lake

abstract type LakeParameterization <: CryoGrid.Parameterization end

struct SimpleLakeScheme <: LakeParameterization end

Base.@kwdef struct Lake{Tpara<:LakeParameterization,Tsp,Tproc,Tprop} <: CryoGrid.SubSurface{Tproc}
    para::Tpara = SimpleLakeScheme()
    prop::Tprop = CryoGrid.ThermalProperties()
    sp::Tsp = nothing
    proc::Tproc
end

Lake(proc::Tproc; kwargs...) where {Tproc} = Lake(;proc, kwargs...)

CryoGrid.variables(lake::Lake, heat::HeatBalance) = (
    CryoGrid.variables(heat)...,
    Diagnostic(:ρ_w, Scalar, u"kg*m^-3", domain=0..Inf, desc = "density of water with temperature"),
    Diagnostic(:T_ub, Scalar, u"°C"),
)

function CryoGrid.initialcondition!(lake::Lake, heat::HeatBalance, state)
    L = heat.prop.L
    # initialize liquid water content based on temperature
    @inbounds for i in 1:length(state.T)
        θwi = Hydrology.watercontent(lake, state, i)
        state.θw[i] = ifelse(state.T[i] > 0.0, θwi, 0.0)
        state.C[i] = heatcapacity(lake, heat, volumetricfractions(lake, state, i)...)
        state.H[i] = enthalpy(state.T[i], state.C[i], L, state.θw[i])
    end
end

function CryoGrid.diagnosticstep!(sub::Lake, heat::HeatBalance, state)
    Heat.resetfluxes!(sub, heat, state)
    # Evaluate freeze/thaw processes
    Heat.freezethaw!(sub, heat, state)

    # force unfrozen water to Tair
    if all(state.Θw .>= 1.)
        state.T .= state.T_ub
        state.H .= Heat.enthalpy(state.T, state.C, heat.prop.L, state.Θw)
    end

    # Update thermal conductivity
    Heat.thermalconductivity!(sub, heat, state)
    return nothing
end


function CryoGrid.prognosticstep!(::Lake, ::HeatBalance{<:FreezeCurve,<:Heat.Enthalpy}, state)
    Δk = Δ(state.grids.k) # cell sizes
    ΔT = Δ(state.grids.T) # midpoint distances
    # compute internal fluxes and non-linear diffusion assuming boundary fluxes have been set
    Numerics.nonlineardiffusion!(state.∂H∂t, state.jH, state.T, ΔT, state.k, Δk)
    return nothing
end

function CryoGrid.interact!(
    top::Top,
    bc::HeatBC,
    lake::Lake,
    heat::HeatBalance,
    stop,
    slake
)
    slake.T_ub = CryoGrid.boundaryvalue(bc, top, heat, lake, stop, slake)
    # boundary flux
    slake.jH[1] += CryoGrid.boundaryflux(bc, top, heat, lake, stop, slake)
    return nothing
end


end