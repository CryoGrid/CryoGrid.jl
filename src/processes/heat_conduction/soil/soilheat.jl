""" Defined variables for heat conduction (enthalpy) on soil layer. """
variables(soil::Soil, heat::Heat{UT"J"}) = (
    Prognostic(:H, Float"J", OnGrid(Cells)),
    Diagnostic(:T, Float"K", OnGrid(Cells)),
    Diagnostic(:C, Float"J/(K*m^3)", OnGrid(Cells)),
    Diagnostic(:k, Float"W/(m*K)", OnGrid(Edges)),
    Diagnostic(:kc, Float"W/(m*K)", OnGrid(Cells)),
)

function initialcondition!(soil::Soil, heat::Heat{UT"J"}, state)
    interpolateprofile!(heat.profile, state)
    L = ρ(heat)*Lsl(heat)
    @. state.C = heatCapacity(soil.hcparams, state.θw, state.θl, state.θm, state.θo)
    @. state.H = enthalpy(heat, state.T, state.C, state.θl, L)
    # k lies on the boundary grid, so we have to regrid the soil properties
    @. state.kc = thermalConductivity(soil.tcparams, state.θw, state.θl, state.θm, state.θo)
    regrid!(state.k, state.kc, state.grids.kc, state.grids.k, Linear(), Flat())
end

function prognosticstep!(soil::Soil, heat::Heat{UT"J"}, state)
    Δk = Δ(state.grids.k) # cell sizes
    ΔT = Δ(state.grids.T)
    heatconduction!(state.T,ΔT,state.k,Δk,state.dH)
end

include("soilheat_freewater.jl")
