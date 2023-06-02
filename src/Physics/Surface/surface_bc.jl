Base.@kwdef struct SurfaceEnergyWaterFluxes{TH<:HeatBC,TW<:WaterBC} <: BoundaryProcess{Union{WaterBalance,HeatBalance}}
    heat::TH
    water::TW
end

function CryoGrid.updatestate!(top::Top, bc::SurfaceEnergyWaterFluxes, state)
    updatestate!(top, bc.water, state)
    updatestate!(top, bc.heat, state)
end

function CryoGrid.interact!(top::Top, bc::SurfaceEnergyWaterFluxes, sub::SubSurface, proc::SubSurfaceProcess, stop, ssub)
    interact!(top, bc.water, sub, proc, stop, ssub)
    interact!(top, bc.heat, sub, proc, stop, ssub)
end

function CryoGrid.computefluxes!(top::Top, bc::SurfaceEnergyWaterFluxes, state)
    computefluxes!(top, bc.water, state)
    computefluxes!(top, bc.heat, state)
end
