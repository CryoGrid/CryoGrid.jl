# Water/heat coupling

const WaterHeatBC{TWater,THeat} = Coupled2{TWater,THeat} where {TWater<:WaterBC,THeat<:HeatBC}

WaterHeatBC(waterbc::WaterBC, heatbc::HeatBC) = Coupled(waterbc, heatbc)

"""
    advectiveflux(jw, T₁, T₂, cw, L)

Computes the advective energy flux between grid cells with temperatures `T₁` and `T₂`
given the heat capacity of water `cw` and latent heat of fusion `L`.
"""
advectiveflux(jw, T₁, T₂, cw, L) = jw*(cw*T₁*(jw > zero(jw)) + cw*T₂*(jw < zero(jw)) + L)

"""
    water_energy_advection!(::SubSurface, ::Coupled(WaterBalance, HeatBalance), state)

Adds advective energy fluxes for all internal grid cell faces.
"""
function water_energy_advection!(sub::SubSurface, ps::Coupled(WaterBalance, HeatBalance), state)
    water, heat = ps
    cw = heatcapacity_water(sub)
    @inbounds for i in 2:length(state.jw)-1
        let jw = state.jw[i],
            T₁ = state.T[i-1],
            T₂ = state.T[i],
            L = heat.prop.L;
            jH_w = advectiveflux(jw, T₁, T₂, cw, L)
            state.jH[i] += jH_w
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
    if heat.advection
        T_ub = getscalar(stop.T_ub)
        Ts = ssub.T[1]
        cw = heatcapacity_water(sub)
        L = heat.prop.L
        jw = ssub.jw[1]
        jH_w = advectiveflux(jw, T_ub, Ts, cw, L)
        ssub.jH[1] += jH_w
    end
    return nothing
end
function CryoGrid.interact!(sub::SubSurface, heat::HeatBalance, bot::Bottom, bc::WaterBC, ssub, sbot)
    if heat.advection
        Ts = ssub.T[end]
        T_lb = getscalar(sbot.T_ub)
        cw = heatcapacity_water(sub)
        L = heat.prop.L
        jw = ssub.jw[end]
        jH_w = advectiveflux(jw, Ts, T_lb, cw, L)
        ssub.jH[end] += jH_w
    end
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
    if p1[2].advection && p2[2].advection
        # advective energy flux
        T₁ = state1.T[end]
        T₂ = state2.T[1]
        jw = state1.jw[end]
        L = p1[2].prop.L
        cw1 = heatcapacity_water(sub1)
        cw2 = heatcapacity_water(sub2)
        cw = (jw >= zero(jw))*cw1 + (jw < zero(jw))*cw2
        jH_w = advectiveflux(jw, T₁, T₂, cw, L)
        state1.jH[end] = state2.jH[1] += jH_w
    end
end

# Flux calculation
function CryoGrid.computefluxes!(sub::SubSurface, ps::Coupled(WaterBalance, HeatBalance), state)
    water, heat = ps
    CryoGrid.computefluxes!(sub, water, state)
    if heat.advection
        water_energy_advection!(sub, ps, state)
    end
    CryoGrid.computefluxes!(sub, heat, state)
end

function heatcapacity_water(sub::SubSurface)
    @unpack ch_w = thermalproperties(sub)
    return ch_w
end
