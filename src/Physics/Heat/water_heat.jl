# Water/heat coupling
CryoGrid.initialcondition!(sub::SubSurface, ps::Coupled(WaterBalance, HeatBalance), state) = CryoGrid.diagnosticstep!(sub, ps, state)
function CryoGrid.diagnosticstep!(sub::SubSurface, ps::Coupled(WaterBalance, HeatBalance), state)
    water, heat = ps
    # Reset fluxes
    Hydrology.resetfluxes!(sub, water, state)
    # Compute water diagnostics
    Hydrology.watercontent!(sub, water, state)
    # HeatBalance diagnostics
    Heat.diagnosticstep!(sub, heat, state)
    # then hydraulic conductivity (requires liquid water content from heat conduction)
    Hydrology.hydraulicconductivity!(sub, water, state)
end
function CryoGrid.prognosticstep!(sub::SubSurface, ps::Coupled(WaterBalance, HeatBalance), state)
    water, heat = ps
    CryoGrid.prognosticstep!(sub, water, state)
    # heat flux due to change in water content
    @. state.∂H∂t += state.∂θwi∂t*(state.T*(heat.prop.cw - heat.prop.ci) + heat.prop.L)
    CryoGrid.prognosticstep!(sub, heat, state)
end
