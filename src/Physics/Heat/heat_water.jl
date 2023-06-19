# Water/heat coupling

const WaterHeatBC{TWater,THeat} = Coupled2{TWater,THeat} where {TWater<:WaterBC,THeat<:HeatBC}

WaterHeatBC(waterbc::WaterBC, heatbc::HeatBC) = Coupled(waterbc, heatbc)

"""
    advectiveflux(jw, C, T)

Advective energy flux: `jw*C*T` (TODO: should latent heat also be accounted for?)
"""
advectiveflux(jw, C, T) = jw*C*T

"""
    energyadvection!(::SubSurface, ::Coupled(WaterBalance, HeatBalance), state)

Adds advective energy fluxes for all internal grid cell faces.
"""
function energyadvection!(::SubSurface, ::Coupled(WaterBalance, HeatBalance), state)
    @inbounds for i in 2:length(state.jw)-1
        let jw = state.jw[i],
            idx = ifelse(jw < 0, i, i-1),
            C = state.C[idx],
            T = state.T[idx];
            state.jH[i] += advectiveflux(jw, C, T)
        end
    end
end

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

# Flux calculation
function CryoGrid.computefluxes!(sub::SubSurface, ps::Coupled(WaterBalance, HeatBalance), state)
    water, heat = ps
    CryoGrid.computefluxes!(sub, water, state)
    energyadvection!(sub, ps, state)
    CryoGrid.computefluxes!(sub, heat, state)
end
