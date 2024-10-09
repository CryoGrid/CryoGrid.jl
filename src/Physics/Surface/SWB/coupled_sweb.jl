"""
Type alias for `WaterHeatBC{TSWB,TSEB} where {TSWB<:SurfaceWaterBalance,TSEB<:SurfaceEnergyBalance}`.
"""
const SurfaceWaterEnergyBalance{TSWB,TSEB} = WaterHeatBC{TSWB,TSEB} where {TSWB<:SurfaceWaterBalance,TSEB<:SurfaceEnergyBalance}

function ET!(::Top, ::SurfaceWaterBalance, state)
    @setscalar state.dET += state.jw_ET[1]*area(state.grid)
end

CryoGrid.variables(top::Top, bc::SurfaceWaterEnergyBalance) = (
    CryoGrid.variables(top, bc[1])...,
    CryoGrid.variables(top, bc[2])...,
    Prognostic(:ET, Scalar, u"m^3"),
)

function CryoGrid.computeprognostic!(top::Top, bc::SurfaceWaterEnergyBalance, state)
    swb, seb = bc
    computeprognostic!(top, swb, state)
    ET!(top, swb, state)
    computeprognostic!(top, seb, state)
end

function CryoGrid.interact!(top::Top, bc::SurfaceWaterEnergyBalance, sub::SubSurface, ps::Coupled(WaterBalance, HeatBalance), stop, ssub)
    swb, seb = bc
    water, heat = ps
    # first compute SEB interaction for heat
    CryoGrid.interact!(top, seb, sub, heat, stop, ssub)
    # then ET interaction
    Hydrology.interact_ET!(top, bc, sub, water, stop, ssub)
    # then SWB interaction for water...
    CryoGrid.interact!(top, swb, sub, water, stop, ssub)
    # ..and for heat
    CryoGrid.interact!(top, swb, sub, heat, stop, ssub)
end

function CryoGrid.interact!(
    top::Top,
    bc::SurfaceWaterEnergyBalance,
    snow::Snowpack,
    ps::Coupled(SnowMassBalance, WaterBalance, HeatBalance),
    stop,
    ssnow
)
    swb, seb = bc
    mass, water, heat = ps
    # snow mass interaction
    CryoGrid.interact!(top, swb, snow, mass, stop, ssnow)
    # water/heat interactions
    CryoGrid.interact!(top, bc, snow, Coupled(water, heat), stop, ssnow)
end

function Hydrology.interact_ET!(
    ::Top,
    ::SurfaceWaterEnergyBalance,
    sub::SubSurface,
    water::WaterBalance{<:BucketScheme,<:DampedET},
    stop,
    ssub
)
    # propagate surface latent heat flux to next layer
    ssub.Qe .= stop.Qe
    # compute ET fluxes for subsurface layer
    Hydrology.evapotranspirative_fluxes!(sub, water, ssub)
end
