"""
Type alias for the implicit enthalpy formulation of HeatBalance.
"""
const HeatBalanceImplicit = HeatBalance{<:EnthalpyImplicit}

# CryoGrid methods

CryoGrid.variables(heat::HeatBalanceImplicit) = (
    Prognostic(:H, OnGrid(Cells), u"J/m^3"),
    Diagnostic(:T, OnGrid(Cells), u"°C"),
    # coefficients and cache variables for implicit diffusion operator
    Diagnostic(:DT_an, OnGrid(Cells)),
    Diagnostic(:DT_as, OnGrid(Cells)),
    Diagnostic(:DT_ap, OnGrid(Cells)),
    Diagnostic(:DT_bp, OnGrid(Cells)),
    Heat.heatvariables(heat)...,
)

function CryoGrid.computediagnostic!(
    sub::SubSurface,
    heat::HeatBalanceImplicit,
    state
)
    # Evaluate freeze/thaw processes
    freezethaw!(freezecurve(sub), sub, processes(sub), state)
    # Update thermal conductivity
    thermalconductivity!(sub, state)
    # Compute diffusion coefficients
    an = state.DT_an
    as = state.DT_as
    ap = state.DT_ap
    k = state.k
    dx = Δ(cells(state.grid))
    dxp = Δ(state.grid)
    prefactors!(ap, an, as, k, dx, dxp)
    return nothing
end

function CryoGrid.boundaryflux(::Dirichlet, bc::HeatBC, top::Top, heat::HeatBalanceImplicit, sub::SubSurface, stop, ssub)
    T_ub = boundaryvalue(bc, stop)
    k = ssub.k[1]
    Δk = CryoGrid.thickness(sub, ssub, first)
    return 2*T_ub*k / Δk
end
function CryoGrid.boundaryflux(::Dirichlet, bc::HeatBC, bot::Bottom, heat::HeatBalanceImplicit, sub::SubSurface, sbot, ssub)
    T_lb = boundaryvalue(bc, sbot)
    k = ssub.k[end]
    Δk = CryoGrid.thickness(sub, ssub, last)
    return 2*T_lb*k / Δk
end

function CryoGrid.interact!(top::Top, bc::HeatBC, sub::SubSurface, heat::HeatBalanceImplicit, stop, ssub)
    Δk = CryoGrid.thickness(sub, ssub, first)
    jH_top = boundaryflux(bc, top, heat, sub, stop, ssub)
    k = ssub.k[1]
    ssub.DT_bp[1] += jH_top / Δk
    ssub.DT_ap[1] += apbc(CryoGrid.BCKind(bc), k, Δk)
    return nothing
end
function CryoGrid.interact!(sub::SubSurface, heat::HeatBalanceImplicit, bot::Bottom, bc::HeatBC, ssub, sbot)
    Δk = CryoGrid.thickness(sub, ssub, last)
    jH_bot = boundaryflux(bc, bot, heat, sub, sbot, ssub)
    k = ssub.k[1]
    ssub.DT_bp[end] += jH_bot / Δk
    ssub.DT_ap[end] += apbc(CryoGrid.BCKind(bc), k, Δk)
    return nothing
end
function CryoGrid.interact!(sub1::SubSurface, ::HeatBalanceImplicit, sub2::SubSurface, ::HeatBalanceImplicit, s1, s2)
    Δk₁ = CryoGrid.thickness(sub1, s1, last)
    Δk₂ = CryoGrid.thickness(sub2, s2, first)
    Δz = CryoGrid.midpoint(sub2, s2, first) - CryoGrid.midpoint(sub1, s1, last)
    # thermal conductivity between cells
    k = s1.k[end] = s2.k[1] =
        @inbounds let k₁ = s1.kc[end],
            k₂ = s2.kc[1],
            Δ₁ = Δk₁[end],
            Δ₂ = Δk₂[1];
            harmonicmean(k₁, k₂, Δ₁, Δ₂)
        end
    s1.DT_ap[end] += s1.DT_as[end] = k / Δz / Δk₁
    s2.DT_ap[1] += s2.DT_an[1] = k / Δz / Δk₂
    return nothing
end

# do nothing in computefluxes!
CryoGrid.computefluxes!(::SubSurface, ::HeatBalanceImplicit, state) = nothing

function CryoGrid.resetfluxes!(sub::SubSurface, heat::HeatBalanceImplicit, state)
    @inbounds for i in 1:length(state.H)
        state.∂H∂T[i] = 0.0
        state.∂θw∂T[i] = 0.0
        state.DT_ap[i] = 0.0
        state.DT_bp[i] = 0.0
        state.DT_an[i] = 0.0
        state.DT_as[i] = 0.0
    end
end

# implementations

apbc(::Dirichlet, k, Δk) = k / (Δk^2/2)
apbc(::Dirichlet, k, Δk, Δx) = k / Δx / Δk
apbc(::Neumann, k, Δk) = 0

function prefactors!(ap, an, as, k, dx, dxp)
    # loop over grid cells
    @inbounds for i in eachindex(dxp)
        if i == 1
            as[1] = k[2] / dx[1] / dxp[1]
        elseif i == length(dxp)
            an[end] = k[end-1] / dx[end] / dxp[end]
        else
            an[i] = k[i] / dx[i-1] / dxp[i]
            as[i] = k[i+1] / dx[i] / dxp[i]
        end
        ap[i] = an[i] + as[i]
    end
    return nothing
end
