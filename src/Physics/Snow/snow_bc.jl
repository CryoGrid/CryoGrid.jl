# unconditionally apply top interactions for snow
CryoGrid.caninteract(::Top, ::WaterBC, ::Snowpack, ::SnowMassBalance, s1, s2) = true
CryoGrid.caninteract(::Top, ::HeatBC, snow::Snowpack, ::HeatBalance, s1, s2) = true
# unconditionally apply subsurface interactions for water
CryoGrid.caninteract(::Snowpack, ::WaterBalance, ::SubSurface, ::WaterBalance, s1, s2) = true
# only apply heat interactions when snowpack is active
CryoGrid.caninteract(snow::Snowpack, ::HeatBalance, ::SubSurface, ::HeatBalance, s1, s2) = isactive(snow, s2)

# default top interact! for dynamic snow mass balance
function CryoGrid.interact!(
    top::Top,
    sbc::SnowBC,
    snow::Snowpack,
    mass::SnowMassBalance,
    stop,
    ssnow
)
    Snow.accumulation!(top, sbc, snow, mass, stop, ssnow)
    Snow.ablation!(top, sbc, snow, mass, stop, ssnow)
    return nothing
end
# default interact! for heat
function CryoGrid.interact!(
    top::Top,
    bc::HeatBC,
    snow::Snowpack,
    heat::HeatBalance,
    stop,
    ssnow
)
    T_ub = getscalar(stop.T_ub)
    @setscalar ssnow.T_ub = T_ub
    if isactive(snow, ssnow)
        # boundary flux
        ssnow.jH[1] += CryoGrid.boundaryflux(bc, top, heat, snow, stop, ssnow)
    else
        # set snow grid cells temperature to upper boundary temperature
        ssnow.T .= T_ub
    end
    return nothing
end
# default interact! for coupled snow mass, water, and heat
function CryoGrid.interact!(top::Top, bc::WaterHeatBC, snow::Snowpack, ps::CoupledSnowWaterHeat, stop, ssub)
    waterbc, heatbc = bc
    snowmass, water, heat = ps
    interactmaybe!(top, waterbc, snow, snowmass, stop, ssub)
    interactmaybe!(top, waterbc, snow, water, stop, ssub)
    interactmaybe!(top, heatbc, snow, heat, stop, ssub)
end
