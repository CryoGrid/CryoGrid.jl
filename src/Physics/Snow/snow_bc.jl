const SnowBC = BoundaryProcess{T} where {SnowMassBalance<:T<:SubSurfaceProcess}

# always allow Top interactions with Snowpack
CryoGrid.caninteract(::Top, ::WaterBC, ::Snowpack, ::SnowMassBalance, s1, s2) = true
CryoGrid.caninteract(::Top, ::HeatBC, snow::Snowpack, ::HeatBalance, s1, s2) = isactive(snow, s2)

# default top interact! for snow mass balance
function CryoGrid.interact!(
    top::Top,
    sbc::SnowBC,
    snow::Snowpack,
    mass::SnowMassBalance,
    stop,
    ssnow
)
    Snow.accumulate!(top, sbc, snow, mass, stop, ssnow)
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
    @setscalar ssnow.T_ub = getscalar(stop.T_ub)
    # boundary flux
    ssnow.jH[1] += CryoGrid.boundaryflux(bc, top, heat, snow, stop, ssnow)
    return nothing
end
# default interact! for coupled water/heat
function CryoGrid.interact!(top::Top, bc::WaterHeatBC, snow::Snowpack, ps::Coupled(SnowMassBalance, HeatBalance), stop, ssub)
    waterbc, heatbc = bc
    snowmass, heat = ps
    interactmaybe!(top, waterbc, snow, snowmass, stop, ssub)
    interactmaybe!(top, heatbc, snow, heat, stop, ssub)
end
function CryoGrid.interact!(top::Top, bc::WaterHeatBC, snow::Snowpack, ps::Coupled(SnowMassBalance, WaterBalance, HeatBalance), stop, ssub)
    waterbc, heatbc = bc
    snowmass, water, heat = ps
    interactmaybe!(top, waterbc, snow, snowmass, stop, ssub)
    interactmaybe!(top, waterbc, snow, water, stop, ssub)
    interactmaybe!(top, heatbc, snow, heat, stop, ssub)
end
