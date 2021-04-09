""" Defined variables for heat conduction (enthalpy) on soil layer. """
variables(soil::Soil, heat::Heat{u"J"}) = (
    Prognostic(:H, Float"J", OnGrid(Cells)),
    Diagnostic(:T, Float"K", OnGrid(Cells)),
    Diagnostic(:C, Float"J/(K*m^3)", OnGrid(Cells)),
    Diagnostic(:k, Float"W/(m*K)", OnGrid(Edges)),
    Diagnostic(:kc, Float"W/(m*K)", OnGrid(Cells)),
    variables(freezecurve(heat))...,
)

function initialcondition!(soil::Soil, heat::Heat{u"J"}, state)
    interpolateprofile!(heat.profile, state)
    L = heat.params.L
    @. state.C = heatcapacity(soil.hcparams, state.θw, state.θl, state.θm, state.θo)
    @. state.H = enthalpy(state.T, state.C, L, state.θl)
    # k lies on the boundary grid, so we have to regrid the soil properties
    @. state.kc = thermalconductivity(soil.tcparams, state.θw, state.θl, state.θm, state.θo)
    regrid!(state.k, state.kc, state.grids.kc, state.grids.k, Linear(), Flat())
end

function diagnosticstep!(soil::Soil, heat::Heat{u"J"}, state)
    fc! = freezecurve(heat);
    @. state.C = heatcapacity(soil.hcparams, state.θw, state.θl, state.θm, state.θo)
    fc!(soil,heat,state)
    @. state.kc = thermalconductivity(soil.tcparams, state.θw, state.θl, state.θm, state.θo)
    regrid!(state.k, state.kc, state.grids.kc, state.grids.k, Linear(), Flat())
    return nothing # ensure no allocation
end

function prognosticstep!(soil::Soil, heat::Heat{u"J"}, state)
    Δk = Δ(state.grids.k) # cell sizes
    ΔT = Δ(state.grids.T)
    heatconduction!(state.T,ΔT,state.k,Δk,state.dH)
end

include("sfcc.jl")

# Note for future use: harmonic mean of thermal conductivities
# for i=2:N-1
#     kn(i,1) = (dxp(i,1)/(2*dxn(i))*kp(i,1).^-1 + dxp(i-1,1)/(2*dxn(i))*kp(i-1).^-1).^-1;
#     ks(i,1) = (dxp(i,1)/(2*dxs(i))*kp(i,1).^-1 + dxp(i+1,1)/(2*dxs(i))*kp(i+1).^-1).^-1;
# end
