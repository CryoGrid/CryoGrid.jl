Base.@kwdef struct WaterHeatBC{TW<:WaterBC,TH<:HeatBC} <: BoundaryProcess{Union{WaterBalance,HeatBalance}}
    water::TW
    heat::TH
end

CryoGrid.variables(top::Top, bc::WaterHeatBC) = (
    variables(top, bc.water)...,
    variables(top, bc.heat)...,
)

function CryoGrid.updatestate!(top::Top, bc::WaterHeatBC, state)
    updatestate!(top, bc.water, state)
    updatestate!(top, bc.heat, state)
end

function CryoGrid.computefluxes!(top::Top, bc::WaterHeatBC, state)
    computefluxes!(top, bc.water, state)
    computefluxes!(top, bc.heat, state)
end

function CryoGrid.interact!(top::Top, bc::WaterHeatBC, sub::SubSurface, ps::Coupled(WaterBalance, HeatBalance), stop, ssub)
    water, heat = ps
    interactmaybe!(top, bc.water, sub, water, stop, ssub)
    interactmaybe!(top, bc.heat, sub, heat, stop, ssub)
end

function CryoGrid.interact!(top::Top, bc::WaterHeatBC, sub::SubSurface, heat::HeatBalance, stop, ssub)
    interactmaybe!(top, bc.heat, sub, heat, stop, ssub)
end

function CryoGrid.interact!(top::Top, bc::WaterHeatBC, sub::SubSurface, water::WaterBalance, stop, ssub)
    interactmaybe!(top, bc.water, sub, water, stop, ssub)
end

function CryoGrid.caninteract(top::Top, bc::WaterHeatBC, sub::SubSurface, ps::Coupled(WaterBalance, HeatBalance), stop, ssub)
    water, heat = ps
    return caninteract(top, bc.water, sub, water, stop, ssub) || caninteract(top, bc.heat, sub, heat, stop, ssub)
end

function CryoGrid.caninteract(top::Top, bc::WaterHeatBC, sub::SubSurface, heat::HeatBalance, stop, ssub)
    return caninteract(top, bc.heat, sub, heat, stop, ssub)
end

function CryoGrid.caninteract(top::Top, bc::WaterHeatBC, sub::SubSurface, water::WaterBalance, stop, ssub)
    return caninteract(top, bc.water, sub, water, stop, ssub)
end

# SubSurface

function CryoGrid.initialcondition!(sub::SubSurface, ps::Coupled(WaterBalance, HeatBalance), state)
    water, heat = ps
    CryoGrid.initialcondition!(sub, water, state)
    CryoGrid.initialcondition!(sub, heat, state)
end

function CryoGrid.interact!(
    sub1::SubSurface,
    p1::Coupled(WaterBalance, HeatBalance),
    sub2::SubSurface,
    p2::Coupled(WaterBalance, HeatBalance),
    state1,
    state2,
)
    # water interaction
    interact!(sub1, p1[1], sub2, p2[1], state1, state2)
    # heat interaction
    interact!(sub1, p1[2], sub2, p2[2], state1, state2)
end
