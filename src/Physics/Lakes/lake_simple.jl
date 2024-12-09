abstract type LakeParameterization <: CryoGrid.Parameterization end

Base.@kwdef struct SimpleLakeScheme{Thp} <: LakeParameterization
    heat::Thp = ThermalProperties()
end

Base.@kwdef struct Lake{Tpara<:LakeParameterization,Theat<:HeatBalance,Taux} <: CryoGrid.SubSurface
    para::Tpara = SimpleLakeScheme()
    heat::Theat = HeatBalance()
    aux::Taux = nothing
end

function get_lake_ub_idx(T_ub, θw)
    ubc_idx = 1
    if T_ub >= zero(T_ub)
        @inbounds for i in eachindex(θw)
            ubc_idx = i
            if θw[i] < one(eltype(θw))
                return ubc_idx
            end
        end
    end
    return ubc_idx+1
end

# Material properties
Heat.thermalproperties(lake::Lake) = lake.para.heat

Hydrology.watercontent(::Lake, state) = 1.0

CryoGrid.processes(lake::Lake) = lake.heat

CryoGrid.volumetricfractions(::Lake, state, i) = (state.θw[i], 1 - state.θw[i], 0.0)

CryoGrid.variables(::Lake, heat::HeatBalance) = (
    CryoGrid.variables(heat)...,
    Diagnostic(:ρ_w, Scalar, u"kg*m^-3", domain=0..Inf, desc = "density of water with temperature"),
    Diagnostic(:T_ub, Scalar, u"°C"),
    Diagnostic(:ubc_idx, Scalar, NoUnits, Int),
)

function CryoGrid.diagnosticstep!(::Lake, state)
    T_ub = getscalar(state.T_ub)
    @setscalar state.ubc_idx = get_lake_ub_idx(T_ub, state.θw)
    return false
end

function CryoGrid.initialcondition!(lake::Lake, state)
    initialcondition!(lake, lake.heat, state)
    diagnosticstep!(lake, state)
end

function CryoGrid.interact!(top::Top, bc::HeatBC, lake::Lake, heat::HeatBalanceImplicit, stop, slake)
    T_ub = slake.T_ub[1] = getscalar(stop.T_ub)
    ubc_idx = Int(getscalar(slake.ubc_idx))
    # get variables
    an = slake.DT_an
    as = slake.DT_as
    ap = slake.DT_ap
    bp = slake.DT_bp
    k = slake.k
    dx = Δ(cells(slake.grid))
    dxp = Δ(slake.grid)
    # outer lake cell
    bp[1] = T_ub*k[1] / (dxp[1]/2) / dxp[1]
    ap[1] = k[1] / (dxp[1]/2) / dxp[1]
    if ubc_idx > 1
        as[1] = an[1] = zero(eltype(as))
    end
    # inner lake cells
    @inbounds for i in 2:ubc_idx-1
        bp[i] = T_ub*k[i] / dx[i-1] / dxp[i]
        ap[i] = an[i]
        as[i] = zero(eltype(as))
        an[i] = zero(eltype(an))
    end
    return nothing
end
function CryoGrid.interact!(lake::Lake, ::HeatBalanceImplicit, sub::SubSurface, ::HeatBalanceImplicit, slake, ssub)
    Δk₁ = CryoGrid.thickness(lake, slake, last)
    Δk₂ = CryoGrid.thickness(sub, ssub, first)
    Δz = CryoGrid.midpoint(sub, ssub, first) - CryoGrid.midpoint(lake, slake, last)
    ubc_idx = Int(getscalar(slake.ubc_idx))
    # thermal conductivity between cells
    k = slake.k[end] = ssub.k[1] =
        @inbounds let k₁ = slake.kc[end],
            k₂ = ssub.kc[1],
            Δ₁ = Δk₁[end],
            Δ₂ = Δk₂[1];
            harmonicmean(k₁, k₂, Δ₁, Δ₂)
        end
    slake.DT_ap[end] += slake.DT_as[end] = (ubc_idx <= length(slake.θw))*k / Δz / Δk₁
    ssub.DT_ap[1] += ssub.DT_an[1] = k / Δz / Δk₂
    return nothing
end

function CryoGrid.computediagnostic!(sub::Lake, heat::HeatBalance{<:Heat.Diffusion1D{:H}}, state)
    Heat.resetfluxes!(sub, heat, state)
    # Evaluate freeze/thaw processes
    Heat.freezethaw!(freezecurve(sub), sub, processes(sub), state)

    # force unfrozen water to Tair
    if all(state.θw .>= 1.)
        state.T .= state.T_ub
        state.H .= Heat.enthalpy.(state.T, state.C, heat.prop.L, state.θw)
    end

    # Update thermal conductivity
    Heat.thermalconductivity!(sub, state)
    return nothing
end

function CryoGrid.computeprognostic!(::Lake, ::HeatBalance{<:Heat.Diffusion1D{:H}}, state)
    Δk = Δ(state.grids.k) # cell sizes
    ΔT = Δ(state.grids.T) # midpoint distances
    # compute internal fluxes and non-linear diffusion assuming boundary fluxes have been set
    Numerics.nonlineardiffusion!(state.dH, state.jH, state.T, ΔT, state.k, Δk)
    return nothing
end


function CryoGrid.interact!(
    top::Top,
    bc::HeatBC,
    lake::Lake,
    heat::HeatBalance{<:Heat.Diffusion1D{:H}},
    stop,
    slake
)
    @setscalar slake.T_ub = stop.T_ub
    # boundary flux
    slake.jH[1] += CryoGrid.boundaryflux(bc, top, heat, lake, stop, slake)
    return nothing
end