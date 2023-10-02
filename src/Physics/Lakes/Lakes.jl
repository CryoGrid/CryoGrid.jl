module Lakes

using CryoGrid
using CryoGrid.Utils
using CryoGrid.Numerics

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

# Material properties
Heat.thermalproperties(lake::Lake) = lake.prop
Hydrology.watercontent(::Lake, state) = 1.0
CryoGrid.volumetricfractions(::Lake, state, i) = (state.θw[i], 1 - state.θw[i], 0.0)

CryoGrid.variables(lake::Lake, heat::HeatBalance) = (
    CryoGrid.variables(heat)...,
    # Diagnostic(:I_t, OnGrid(Cells), desc="Indicator variable for is thawed (1 or 0)."),
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

function CryoGrid.updatestate!(sub::Lake, heat::HeatBalance, state)
    Heat.resetfluxes!(sub, heat, state)
    # Evaluate freeze/thaw processes
    Heat.freezethaw!(sub, heat, state)

    # force unfrozen water to Tair
    if all(state.θw .>= 1.)
        state.T .= state.T_ub
        state.H .= Heat.enthalpy.(state.T, state.C, heat.prop.L, state.θw)
    end

    # Update thermal conductivity
    Heat.thermalconductivity!(sub, heat, state)
    return nothing
end

function CryoGrid.updatestate!(
    sub::Lake,
    heat::HeatBalanceImplicit,
    state
)
    Heat.resetfluxes!(sub, heat, state)
    # Evaluate freeze/thaw processes
    Heat.freezethaw!(sub, heat, state)
    # Update thermal conductivity
    Heat.thermalconductivity!(sub, heat, state)
    isthawed = true
    @inbounds for i in eachindex(state.θw)
        isthawed = isthawed && state.θw[i] ≈ 1.0
    end
    I_f = 1 - Float64(isthawed)
    # Compute diffusion coefficients
    an = state.DT_an
    as = state.DT_as
    ap = state.DT_ap
    k = state.k
    dx = Δ(cells(state.grid))
    dxp = Δ(state.grid)
    k_inner = @view k[2:end-1]
    dxpn = @view dxp[1:end-1]
    dxps = @view dxp[2:end]
    @. an[2:end] = k_inner / dx / dxpn
    @. as[1:end-1] = (k_inner / dx / dxps)*I_f
    @. ap[1:end-1] += as[1:end-1]
    @. ap[2:end] += an[2:end]
    @. an[2:end] *= I_f
    return nothing
end

CryoGrid.computefluxes!(::Lake, ::HeatBalanceImplicit, state) = nothing

function CryoGrid.computefluxes!(::Lake, ::HeatBalance{<:FreezeCurve,<:Heat.Enthalpy}, state)
    Δk = Δ(state.grids.k) # cell sizes
    ΔT = Δ(state.grids.T) # midpoint distances
    # compute internal fluxes and non-linear diffusion assuming boundary fluxes have been set
    Numerics.nonlineardiffusion!(state.∂H∂t, state.jH, state.T, ΔT, state.k, Δk)
    return nothing
end

function CryoGrid.interact!(top::Top, bc::HeatBC, lake::Lake, heat::HeatBalanceImplicit, stop, slake)
    Δk = CryoGrid.thickness(lake, slake, first)
    jH_top = boundaryflux(bc, top, heat, lake, stop, slake)
    k = slake.k[1]
    slake.DT_bp[1] += jH_top / Δk
    slake.DT_ap[1] += Heat._ap(CryoGrid.BoundaryStyle(bc), k, Δk)
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
    @setscalar slake.T_ub = CryoGrid.boundaryvalue(bc, top, heat, lake, stop, slake)
    # boundary flux
    slake.jH[1] += CryoGrid.boundaryflux(bc, top, heat, lake, stop, slake)
    return nothing
end


end