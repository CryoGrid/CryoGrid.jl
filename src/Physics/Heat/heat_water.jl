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
    water_energy_advection!(::SubSurface, ::Coupled(WaterBalance, HeatBalance), state)

Adds advective energy fluxes for all internal grid cell faces.
"""
function water_energy_advection!(sub::SubSurface, ::Coupled(WaterBalance, HeatBalance), state)
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

# Water/heat interactions
function CryoGrid.interact!(top::Top, bc::WaterBC, sub::SubSurface, heat::HeatBalance, stop, ssub)
    T_ub = getscalar(stop.T_ub)
    Ts = ssub.T[1]
    cw = thermalproperties(sub).ch_w
    jw = ssub.jw[1]
    ssub.jH[1] += advectiveflux(jw, cw, T_ub, Ts)
    return nothing
end
function CryoGrid.interact!(sub::SubSurface, heat::HeatBalance, bot::Bottom, bc::WaterBC, ssub, sbot)
    Ts = ssub.T[end]
    T_ub = getscalar(sbot.T_ub)
    cw = thermalproperties(sub).ch_w
    jw = ssub.jw[end]
    ssub.jH[end] += advectiveflux(jw, cw, Ts, T_ub)
    return nothing
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
    # advective energy flux
    T1 = state1.T[end]
    T2 = state2.T[1]
    jw = state1.jw[end]
    cw1 = thermalproperties(sub1).ch_w
    cw2 = thermalproperties(sub2).ch_w
    cw = (jw > 0)*cw1 + (jw < 0)*cw2
    state1.jH[end] += advectiveflux(jw, cw, T1, T2)
end

# Flux calculation
function CryoGrid.computefluxes!(sub::SubSurface, ps::Coupled(WaterBalance, HeatBalance), state)
    water, heat = ps
    CryoGrid.computefluxes!(sub, water, state)
    CryoGrid.computefluxes!(sub, heat, state)
    water_energy_advection!(sub, ps, state)
end
