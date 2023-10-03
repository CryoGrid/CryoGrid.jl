module Lakes

using CryoGrid
using CryoGrid.Utils
using CryoGrid.Numerics

export Lake

abstract type LakeParameterization <: CryoGrid.Parameterization end

Base.@kwdef struct SimpleLakeScheme{Thp} <: LakeParameterization
    heat::Thp = ThermalProperties()
end

Base.@kwdef struct Lake{Tpara<:LakeParameterization,Theat<:HeatBalance,Taux} <: CryoGrid.SubSurface
    para::Tpara = SimpleLakeScheme()
    heat::Theat = HeatBalance()
    aux::Taux = nothing
end

function get_upper_boundary_index(T_ub, θw)
    ubc_idx = 1
    if T_ub >= zero(T_ub)
        @inbounds for i in eachindex(θw)
            ubc_idx = i
            if θw[i] < one(eltype(θw))
                break
            end
        end
    end
    return ubc_idx
end

# Material properties
Heat.thermalproperties(lake::Lake) = lake.para.heat

Hydrology.watercontent(::Lake, state) = 1.0

CryoGrid.processes(lake::Lake) = lake.heat

CryoGrid.volumetricfractions(::Lake, state, i) = (state.θw[i], 1 - state.θw[i], 0.0)

CryoGrid.variables(lake::Lake, heat::HeatBalance) = (
    CryoGrid.variables(heat)...,
    # Diagnostic(:I_t, OnGrid(Cells), desc="Indicator variable for is thawed (1 or 0)."),
    Diagnostic(:ρ_w, Scalar, u"kg*m^-3", domain=0..Inf, desc = "density of water with temperature"),
    Diagnostic(:T_ub, Scalar, u"°C"),
    Diagnostic(:ubc_idx, Scalar, NoUnits, Int),
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

function CryoGrid.interact!(top::Top, bc::HeatBC, lake::Lake, heat::HeatBalanceImplicit, stop, slake)
    jH_top = boundaryflux(bc, top, heat, lake, stop, slake)
    T_ub = slake.T_ub[1] = getscalar(stop.T_ub)
    ubc_idx = get_upper_boundary_index(T_ub, slake.θw)
    # get variables
    an = slake.DT_an
    as = slake.DT_as
    ap = slake.DT_ap
    bp = slake.DT_bp
    k = slake.k
    dx = Δ(cells(slake.grid))
    dxp = Δ(slake.grid)
    # first lake cell
    bp[1] += jH_top / dxp[1]
    ap[1] += Heat.apbc(CryoGrid.BCKind(bc), k[1], dxp[1])
    # deeper lake cells
    @inbounds for i in 2:ubc_idx
        bp[i] += jH_top / dxp[i]
        ap[i] -= as[i]
        as[i] = zero(eltype(as))
        an[i] = zero(eltype(an))
        ap[i] += Heat.apbc(CryoGrid.BCKind(bc), k[i], dxp[i], dx[i-1])
    end
    return nothing
end
function CryoGrid.interact!(lake::Lake, ::HeatBalanceImplicit, sub::SubSurface, ::HeatBalanceImplicit, slake, ssub)
    Δk₁ = CryoGrid.thickness(lake, slake, last)
    Δk₂ = CryoGrid.thickness(sub, ssub, first)
    Δz = CryoGrid.midpoint(sub, ssub, first) - CryoGrid.midpoint(lake, slake, last)
    # thermal conductivity between cells
    k = slake.k[end] = ssub.k[1] =
        @inbounds let k₁ = slake.kc[end],
            k₂ = ssub.kc[1],
            Δ₁ = Δk₁[end],
            Δ₂ = Δk₂[1];
            harmonicmean(k₁, k₂, Δ₁, Δ₂)
        end
    ubc_idx = get_upper_boundary_index(slake.T_ub[1], slake.θw)
    slake.DT_ap[end] += slake.DT_as[end] = (ubc_idx < length(slake.θw))*k / Δz / Δk₁
    ssub.DT_ap[1] += ssub.DT_an[1] = k / Δz / Δk₂
    return nothing
end

function CryoGrid.updatestate!(sub::Lake, heat::HeatBalance{FreeWater,<:Heat.MOLEnthalpy}, state)
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

function CryoGrid.computefluxes!(::Lake, ::HeatBalance{FreeWater,<:Heat.MOLEnthalpy}, state)
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
    heat::HeatBalance{FreeWater,<:Heat.MOLEnthalpy},
    stop,
    slake
)
    @setscalar slake.T_ub = stop.T_ub
    # boundary flux
    slake.jH[1] += CryoGrid.boundaryflux(bc, top, heat, lake, stop, slake)
    return nothing
end

end