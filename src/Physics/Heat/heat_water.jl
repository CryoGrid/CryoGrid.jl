# Water/heat coupling

const WaterHeatBC{TWater,THeat} = Coupled2{TWater,THeat} where {TWater<:WaterBC,THeat<:HeatBC}

WaterHeatBC(waterbc::WaterBC, heatbc::HeatBC) = Coupled(waterbc, heatbc)

"""
    advectiveflux(jw, cw, T₁, T₂)

Computes the advective energy flux between two grid cells with temperatures `T₁` and `T₂`.
This function assumes that `jw` is positive downward such that a positive temperature gradient
`T₁ - T₂` would result in a positive (downward) energy flux `jH`.
"""
advectiveflux(jw, cw, T₁, T₂) = jw*cw*(T₁ - T₂)*sign(jw)

"""
    energyadvection!(::SubSurface, ::Coupled(WaterBalance, HeatBalance), state)

Adds advective energy fluxes for all internal grid cell faces.
"""
function energyadvection!(sub::SubSurface, ::Coupled(WaterBalance, HeatBalance), state)
    @unpack ch_w = thermalproperties(sub)
    @inbounds for i in 2:length(state.jw)-1
        let jw = state.jw[i],
            T₁ = state.jw[i-1],
            T₂ = state.jw[i];
            state.jH[i] += advectiveflux(jw, ch_w, T₁, T₂)
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
