# type alias for water balance with ET scheme
const WaterBalanceET{TET} = WaterBalance{TFlow,TET} where {TFlow,TET<:Evapotranspiration}

function CryoGrid.interact!(top::Top, seb::SurfaceEnergyBalance, sub::SubSurface, ps::Coupled(WaterBalanceET, HeatBalance), stop, ssub)
    water, heat = ps
    # first interact! for heat
    interact!(top, seb, sub, heat, stop, ssub)
    # then water
    interact!(top, seb, sub, water, stop, ssub)
end

function CryoGrid.interact!(::Top, ::SurfaceEnergyBalance, ::SubSurface, ::WaterBalanceET, stop, ssub)
    # copy latent heat flux to subsurface layer diagnostic variable
    @. ssub.Qe = stop.Qe
end
