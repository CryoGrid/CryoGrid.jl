function CryoGrid.diagnosticstep!(soil::Soil, ps::Coupled2{<:WaterFluxes,<:Heat}, state)
    water, heat = ps
    # first set water/ice diagnostics
    Hydrology.waterice!(soil, water, state)
    # then evaluate diagnosticstep! for Heat
    CryoGrid.diagnosticstep!(soil, heat, state)
    # then hydraulic conductivity
    Hydrology.hydraulicconductivity!(soil, water, state)
end
