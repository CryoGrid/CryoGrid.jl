"""
Type alias for the implicit enthalpy formulation of Heat.
"""
const ImplicitHeat{Tfc} = Heat{Tfc,EnthalpyImplicit} where {Tfc}

CryoGrid.variables(heat::ImplicitHeat) = (
    Prognostic(:H, OnGrid(Cells), u"J/m^3"),
    Diagnostic(:T, OnGrid(Cells), u"°C"),
    Diagnostic(:∂H∂T, OnGrid(Cells), u"J/K/m^3", domain=0..Inf),
    Diagnostic(:∂θw∂T, OnGrid(Cells), domain=0..Inf),
    Diagnostic(:C, OnGrid(Cells), u"J/K/m^3"),
    Diagnostic(:k, OnGrid(Edges), u"W/m/K"),
    Diagnostic(:kc, OnGrid(Cells), u"W/m/K"),
    Diagnostic(:θw, OnGrid(Cells), domain=0..1),
    # coefficients and cache variables for diffusion operator
    Diagnostic(:DT_an, OnGrid(Cells, n -> n-1)),
    Diagnostic(:DT_as, OnGrid(Cells, n -> n-1)),
    Diagnostic(:DT_ap, OnGrid(Cells)),
    Diagnostic(:DT_bp, OnGrid(Cells)),
)
function CryoGrid.diagnosticstep!(
    sub::SubSurface,
    heat::ImplicitHeat,
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
    bp = state.DT_bp
    ap .= 0.0
    bp .= 0.0
    k = state.k
    dx = Δ(cells(state.grid))
    dxp = Δ(state.grid)
    k_inner = @view k[2:end-1]
    dxpn = @view dxp[1:end-1]
    dxps = @view dxp[2:end]
    @. an = k_inner / dx / dxpn
    @. as = k_inner / dx / dxps
    @. ap[1:end-1] += as
    @. ap[2:end] += an
    return nothing # ensure no allocation
end
function CryoGrid.boundaryflux(::Dirichlet, ::HeatBC, ::Top, heat::Heat, sub::SubSurface, stop, ssub)
    T_ub = boundaryvalue(bc, top, heat, sub, stop, ssub)
    k = ssub.k[1]
    Δk = CryoGrid.thickness(sub, ssub, first)
    return 2*T_ub*k / Δk^2
end
function CryoGrid.boundaryflux(::Dirichlet, ::HeatBC, ::Bottom, heat::Heat, sub::SubSurface, sbot, ssub)
    T_lb = boundaryvalue(bc, bot, heat, sub, sbot, ssub)
    k = ssub.k[end]
    Δk = CryoGrid.thickness(sub, ssub, last)
    return 2*T_lb*k / Δk
end
function CryoGrid.interact!(top::Top, bc::HeatBC, sub::SubSurface, heat::ImplicitHeat, stop, ssub)
    Δk = CryoGrid.thickness(sub, ssub, first)
    jH_top = boundaryflux(bc, top, heat, sub, stop, ssub)
    ssub.DT_bp[1] += jH_top / Δk
    return nothing # ensure no allocation
end
function CryoGrid.interact!(sub::SubSurface, heat::ImplicitHeat, bot::Bottom, bc::HeatBC, ssub, sbot)
    Δk = CryoGrid.thickness(sub, ssub, last)
    jH_bot = boundaryflux(bc, bot, heat, sub, sbot, ssub)
    ssub.DT_bp[end] += jH_bot / Δk
    return nothing # ensure no allocation
end
function CryoGrid.interact!(sub1::SubSurface, ::ImplicitHeat, sub2::SubSurface, ::ImplicitHeat, s1, s2)
    Δk₁ = CryoGrid.thickness(sub1, s1, last)
    Δk₂ = CryoGrid.thickness(sub2, s2, first)
    Δz = CryoGrid.midpoints(sub2, s2, first) - CryoGrid.midpoints(sub1, s1, last)
    # thermal conductivity between cells
    k = s1.k[end] = s2.k[1] =
        @inbounds let k₁ = s1.kc[end],
            k₂ = s2.kc[1],
            Δ₁ = Δk₁[end],
            Δ₂ = Δk₂[1];
            harmonicmean(k₁, k₂, Δ₁, Δ₂)
        end
    s1.ap[end] += s1.as[end] = k / Δz / Δk₁
    s2.ap[1] += s2.an[1] = k / Δz / Δk₂
    return nothing # ensure no allocation
end
# do nothing in prognostic step
CryoGrid.prognosticstep!(::SubSurface, ::ImplicitHeat, state) where {Tfc} = nothing
