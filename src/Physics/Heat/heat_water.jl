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

# interactions
function CryoGrid.interact!(top::Top, bc::WaterHeatBC, sub::SubSurface, ps::Coupled(WaterBalance, HeatBalance), stop, ssub)
    water, heat = ps
    interact!(top, bc.water, sub, water, stop, ssub)
    interact!(top, bc.heat, sub, heat, stop, ssub)
end

function CryoGrid.interact!(top::Top, bc::WaterHeatBC, sub::SubSurface, heat::HeatBalance, stop, ssub)
    interact!(top, bc.heat, sub, heat, stop, ssub)
end

# Water/heat coupling
function CryoGrid.initialcondition!(sub::SubSurface, ps::Coupled(WaterBalance, HeatBalance), state)
    water, heat = ps
    CryoGrid.initialcondition!(sub, water, state)
    CryoGrid.initialcondition!(sub, heat, state)
end

# interact!
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
