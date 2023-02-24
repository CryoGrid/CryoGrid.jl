"""
Type alias for the implicit enthalpy formulation of HeatBalance.
"""
const HeatBalanceImplicit{Tfc} = HeatBalance{Tfc,<:EnthalpyImplicit} where {Tfc<:FreezeCurve}

function Heat.resetfluxes!(sub::SubSurface, heat::HeatBalanceImplicit, state)
    @inbounds for i in 1:length(state.H)
        state.∂H∂T[i] = 0.0
        state.∂θw∂T[i] = 0.0
        state.DT_ap[i] = 0.0
        state.DT_bp[i] = 0.0
        state.DT_an[i] = 0.0
        state.DT_as[i] = 0.0
    end
end
CryoGrid.variables(::HeatBalanceImplicit) = (
    Prognostic(:H, OnGrid(Cells), u"J/m^3"),
    Diagnostic(:T, OnGrid(Cells), u"°C"),
    Diagnostic(:∂H∂T, OnGrid(Cells), u"J/K/m^3", domain=0..Inf),
    Diagnostic(:∂θw∂T, OnGrid(Cells), domain=0..Inf),
    Diagnostic(:C, OnGrid(Cells), u"J/K/m^3"),
    Diagnostic(:k, OnGrid(Edges), u"W/m/K"),
    Diagnostic(:kc, OnGrid(Cells), u"W/m/K"),
    Diagnostic(:θw, OnGrid(Cells), domain=0..1),
    # coefficients and cache variables for diffusion operator
    Diagnostic(:DT_an, OnGrid(Cells)),
    Diagnostic(:DT_as, OnGrid(Cells)),
    Diagnostic(:DT_ap, OnGrid(Cells)),
    Diagnostic(:DT_bp, OnGrid(Cells)),
)
function CryoGrid.diagnosticstep!(
    sub::SubSurface,
    heat::HeatBalanceImplicit,
    state
)
    resetfluxes!(sub, heat, state)
    # Evaluate freeze/thaw processes
    freezethaw!(sub, heat, state)
    # Update thermal conductivity
    thermalconductivity!(sub, heat, state)
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
    @. as[1:end-1] = k_inner / dx / dxps
    @. ap[1:end-1] += as[1:end-1]
    @. ap[2:end] += an[2:end]
    return nothing
end
function CryoGrid.boundaryflux(::Dirichlet, bc::HeatBC, top::Top, heat::HeatBalanceImplicit, sub::SubSurface, stop, ssub)
    T_ub = boundaryvalue(bc, top, heat, sub, stop, ssub)
    k = ssub.k[1]
    Δk = CryoGrid.thickness(sub, ssub, first)
    return 2*T_ub*k / Δk
end
function CryoGrid.boundaryflux(::Dirichlet, bc::HeatBC, bot::Bottom, heat::HeatBalanceImplicit, sub::SubSurface, sbot, ssub)
    T_lb = boundaryvalue(bc, bot, heat, sub, sbot, ssub)
    k = ssub.k[end]
    Δk = CryoGrid.thickness(sub, ssub, last)
    return 2*T_lb*k / Δk
end
_ap(::Dirichlet, k, Δk) = k / (Δk^2/2)
_ap(::Neumann, k, Δk) = 0
function CryoGrid.interact!(top::Top, bc::HeatBC, sub::SubSurface, heat::HeatBalanceImplicit, stop, ssub)
    Δk = CryoGrid.thickness(sub, ssub, first)
    jH_top = boundaryflux(bc, top, heat, sub, stop, ssub)
    k = ssub.k[1]
    ssub.DT_bp[1] += jH_top / Δk
    ssub.DT_ap[1] += _ap(CryoGrid.BoundaryStyle(bc), k, Δk)
    return nothing
end
function CryoGrid.interact!(sub::SubSurface, heat::HeatBalanceImplicit, bot::Bottom, bc::HeatBC, ssub, sbot)
    Δk = CryoGrid.thickness(sub, ssub, last)
    jH_bot = boundaryflux(bc, bot, heat, sub, sbot, ssub)
    k = ssub.k[1]
    ssub.DT_bp[end] += jH_bot / Δk
    ssub.DT_ap[end] += _ap(CryoGrid.BoundaryStyle(bc), k, Δk)
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
# do nothing in prognostic step
CryoGrid.prognosticstep!(::SubSurface, ::HeatBalanceImplicit, state) = nothing
