"""
Type alias for `WaterHeatBC{TSWB,TSEB} where {TSWB<:SurfaceWaterBalance,TSEB<:SurfaceEnergyBalance}`.
"""
const SurfaceWaterEnergyBalance{TSWB,TSEB} = WaterHeatBC{TSWB,TSEB} where {TSWB<:SurfaceWaterBalance,TSEB<:SurfaceEnergyBalance}

CryoGrid.variables(top::Top, bc::SurfaceWaterEnergyBalance) = (
    CryoGrid.variables(top, bc[1])...,
    CryoGrid.variables(top, bc[2])...,
    Prognostic(:ET, Scalar, u"m^3"),
)

function interact!(top::Top, bc::SurfaceWaterEnergyBalance, sub::SubSurface, ps::Coupled(WaterBalance, HeatBalance), stop, ssub)
    swb, seb = bc.proc
    water, heat = ps
    # first compute SEB interaction for heat
    interact!(top, seb, sub, heat, stop, ssub)
    # then ET interaction
    Hydrology.interact_ET!(top, bc, sub, water, stop, ssub)
    # then SWB interaction for water...
    interact!(top, swb, sub, water, stop, ssub)
    # ..and for heat
    interact!(top, swb, sub, heat, stop, ssub)
end

function interact!(
    top::Top,
    bc::SurfaceWaterEnergyBalance,
    snow::Snowpack,
    ps::Coupled(SnowMassBalance, WaterBalance, HeatBalance),
    stop,
    ssnow
)
    swb, seb = bc.proc
    mass, water, heat = ps
    # snow mass interaction
    interact!(top, water, snow, mass, stop, ssub)
    # water/heat interactions
    interact!(top, bc, snow, Coupled(water, heat), stop, ssub)
end
