function heatcapacity(params::SoilParams, totalWater, liquidWater, mineral, organic)
    @unpack cw,co,cm,ca,ci = params.hc
    let air = 1.0 - totalWater - mineral - organic,
        ice = totalWater - liquidWater,
        water = liquidWater;
        water*cw + ice*ci + mineral*cm + organic*co + air*ca
    end
end

function thermalconductivity(params::SoilParams, totalWater, liquidWater, mineral, organic)
    @unpack kw,ko,km,ka,ki = params.tc
    let air = 1.0 - totalWater - mineral - organic,
        ice = totalWater - liquidWater,
        water = liquidWater;
        (water*kw^0.5 + ice*ki^0.5 + mineral*km^0.5 + organic*ko^0.5 + air*ka^0.5)^2
    end
end

""" Heat capacity for soil layer """
function heatcapacity!(soil::Soil, ::Heat, state)
    @. state.C = heatcapacity(soil.params, state.θw, state.θl, state.θm, state.θo)
end

""" Thermal conductivity for soil layer """
function thermalconductivity!(soil::Soil, ::Heat, state)
    @. state.kc = thermalconductivity(soil.params, state.θw, state.θl, state.θm, state.θo)
end

include("sfcc.jl")

""" Variable definitions for heat conduction (enthalpy) on soil layer. """
CryoGrid.variables(soil::Soil, heat::Heat{:H}) = (
    Prognostic(:H, Float"J/m^3", OnGrid(Cells)),
    Diagnostic(:T, Float"K", OnGrid(Cells)),
    Diagnostic(:C, Float"J//K*/m^3", OnGrid(Cells)),
    Diagnostic(:Ceff, Float"J/K/m^3", OnGrid(Cells)),
    Diagnostic(:k, Float"W/m/K", OnGrid(Edges)),
    Diagnostic(:kc, Float"W//m/K", OnGrid(Cells)),
    CryoGrid.variables(soil, heat, freezecurve(heat))...,
)
""" Variable definitions for heat conduction (partitioned enthalpy) on soil layer. """
CryoGrid.variables(soil::Soil, heat::Heat{(:Hₛ,:Hₗ)}) = (
    Prognostic(:Hₛ, Float"J/m^3", OnGrid(Cells)),
    Prognostic(:Hₗ, Float"J/m^3", OnGrid(Cells)),
    Diagnostic(:dH, Float"J/s/m^3", OnGrid(Cells)),
    Diagnostic(:H, Float"J", OnGrid(Cells)),
    Diagnostic(:T, Float"K", OnGrid(Cells)),
    Diagnostic(:C, Float"J/K/m^3", OnGrid(Cells)),
    Diagnostic(:Ceff, Float"J/K/m^3", OnGrid(Cells)),
    Diagnostic(:dθdT, Float"m/m", OnGrid(Cells)),
    Diagnostic(:k, Float"W/m/K", OnGrid(Edges)),
    Diagnostic(:kc, Float"W/m/K", OnGrid(Cells)),
    CryoGrid.variables(soil, heat, freezecurve(heat))...,
)
""" Variable definitions for heat conduction (temperature) on soil layer. """
CryoGrid.variables(soil::Soil, heat::Heat{:T}) = (
    Prognostic(:T, Float"K", OnGrid(Cells)),
    Diagnostic(:H, Float"J/m^3", OnGrid(Cells)),
    Diagnostic(:dH, Float"J/s/m^3", OnGrid(Cells)),
    Diagnostic(:C, Float"J/K/m^3", OnGrid(Cells)),
    Diagnostic(:Ceff, Float"J/K/m^3", OnGrid(Cells)),
    Diagnostic(:k, Float"W/m/K", OnGrid(Edges)),
    Diagnostic(:kc, Float"W/m/K", OnGrid(Cells)),
    CryoGrid.variables(soil, heat, freezecurve(heat))...,
)

""" Initial condition for heat conduction (all state configurations) on soil layer. """
function CryoGrid.initialcondition!(soil::Soil, heat::Heat, state)
    interpolateprofile!(heat.profile, state)
    L = heat.params.L
    @. state.C = heatcapacity(soil.params, state.θw, state.θl, state.θm, state.θo)
    @. state.H = enthalpy(state.T, state.C, L, state.θl)
end

""" Initial condition for heat conduction (all state configurations) on soil layer. """
function CryoGrid.initialcondition!(soil::Soil, heat::Heat{U,<:SFCC}, state) where U
    interpolateprofile!(heat.profile, state)
    L = heat.params.L
    sfcc = freezecurve(heat)
    state.θl .= sfcc.f.(state.T, sfccparams(sfcc.f, soil, heat, state)...)
    @. state.C = heatcapacity(soil.params, state.θw, state.θl, state.θm, state.θo)
    @. state.H = enthalpy(state.T, state.C, L, state.θl)
end

""" Diagonstic step for heat conduction (all state configurations) on soil layer. """
function CryoGrid.initialcondition!(soil::Soil, heat::Heat{(:Hₛ,:Hₗ),<:SFCC}, state)
    interpolateprofile!(heat.profile, state)
    L = heat.params.L
    sfcc = freezecurve(heat)
    state.θl .= sfcc.f.(state.T, sfccparams(sfcc.f, soil, heat, state)...)
    @. state.C = heatcapacity(soil.params, state.θw, state.θl, state.θm, state.θo)
    @. state.Hₛ = (state.T - 273.15)*state.C
    @. state.Hₗ = state.θl*L
end

""" Diagonstic step for heat conduction (all state configurations) on soil layer. """
function CryoGrid.diagnosticstep!(soil::Soil, heat::Heat, state)
    # Reset energy flux to zero; this is redundant when H is the prognostic variable
    # but necessary when it is not.
    @. state.dH = zero(eltype(state.dH))
    # Evaluate the freeze curve (updates T, C, and θl)
    fc! = freezecurve(heat);
    fc!(soil,heat,state)
    # Update thermal conductivity
    @. state.kc = thermalconductivity(soil.params, state.θw, state.θl, state.θm, state.θo)
    # Interpolate thermal conductivity to boundary grid
    regrid!(state.k, state.kc, state.grids.kc, state.grids.k, Linear(), Flat())
    # TODO: harmonic mean of thermal conductivities (in MATLAB code)
    # for i=2:N-1
    #     kn(i,1) = (dxp(i,1)/(2*dxn(i))*kp(i,1).^-1 + dxp(i-1,1)/(2*dxn(i))*kp(i-1).^-1).^-1;
    #     ks(i,1) = (dxp(i,1)/(2*dxs(i))*kp(i,1).^-1 + dxp(i+1,1)/(2*dxs(i))*kp(i+1).^-1).^-1;
    # end
    return nothing # ensure no allocation
end

""" Prognostic step for heat conduction (enthalpy) on soil layer. """
function CryoGrid.prognosticstep!(::Soil, ::Heat{:H}, state)
    Δk = Δ(state.grids.k) # cell sizes
    ΔT = Δ(state.grids.T)
    # Diffusion on non-boundary cells
    heatconduction!(state.dH,state.T,ΔT,state.k,Δk)
end

""" Prognostic step for heat conduction (partitioned enthalpy) on soil layer."""
function CryoGrid.prognosticstep!(::Soil, heat::Heat{(:Hₛ,:Hₗ)}, state)
    Δk = Δ(state.grids.k) # cell sizes
    ΔT = Δ(state.grids.T)
    # Diffusion on non-boundary cells
    heatconduction!(state.dH,state.T,ΔT,state.k,Δk)
    let L = heat.params.L;
        @. state.dHₛ = state.dH / (L/state.C*state.dθdT + 1)
        # This could also be expressed via a mass matrix with 1
        # in the upper right block diagonal. But this is easier.
        @. state.dHₗ = state.dH - state.dHₛ
    end
end

""" Prognostic step for heat conduction (temperature) on soil layer. """
function CryoGrid.prognosticstep!(::Soil, ::Heat{:T}, state)
    Δk = Δ(state.grids.k) # cell sizes
    ΔT = Δ(state.grids.T)
    # Diffusion on non-boundary cells
    heatconduction!(state.dH,state.T,ΔT,state.k,Δk)
    # Compute temperature flux by dividing by C_eff;
    # C_eff should be computed by the freeze curve.
    @inbounds @. state.dT = state.dH / state.Ceff
    return nothing
end
