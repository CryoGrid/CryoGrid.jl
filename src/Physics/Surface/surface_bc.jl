Base.@kwdef struct WaterHeatBC{TW<:WaterBC,TH<:HeatBC} <: BoundaryProcess{Union{WaterBalance,HeatBalance}}
    water::TW
    heat::TH
end

function CryoGrid.updatestate!(top::Top, bc::WaterHeatBC, state)
    updatestate!(top, bc.water, state)
    updatestate!(top, bc.heat, state)
end

function CryoGrid.interact!(top::Top, bc::WaterHeatBC, sub::SubSurface, ps::Coupled(WaterBalance, HeatBalance), stop, ssub)
    interact!(top, bc.water, sub, ps[1], stop, ssub)
    interact!(top, bc.heat, sub, ps[2], stop, ssub)
end

function CryoGrid.computefluxes!(top::Top, bc::WaterHeatBC, state)
    computefluxes!(top, bc.water, state)
    computefluxes!(top, bc.heat, state)
end
