# Boundary condition type aliases
const Dirichlet = CryoGrid.Dirichlet
const Neumann = CryoGrid.Neumann
const WaterBC = BoundaryProcess{T} where {WaterBalance<:T<:SubSurfaceProcess}

# Constant boundary constructors
ConstantInfiltration(value::Quantity) = ConstantBC(WaterBalance, Neumann, uconvert(u"m/s", value))
ConstantInfiltration(value) = ConstantBC(WaterBalance, Neumann, value)
ImpermeableBoundary() = ConstantBC(WaterBalance, Neumann, 0.0u"m/s")

function balancefluxes!(top::Top, bc::WaterBC, sub::SubSurface, water::WaterBalance, stop, ssub)
    θw = ssub.θw[1]
    θwi = ssub.θwi[1]
    θsat = ssub.θsat[1]
    sat = ssub.sat[1]
    Δz = CryoGrid.thickness(sub, ssub, first)
    jw = ssub.jw_v[1] + ssub.jw_ET[1]
    ssub.jw[1] = limit_upper_flux(water, jw*stop.dt, θw, θwi, θsat, sat, Δz)/stop.dt
end

function balancefluxes!(sub::SubSurface, water::WaterBalance, bot::Bottom, bc::WaterBC, ssub, sbot)
    θw = ssub.θw[end]
    θwi = ssub.θwi[end]
    θsat = ssub.θsat[end]
    sat = ssub.sat[end]
    Δz = CryoGrid.thickness(sub, ssub, last)
    jw = ssub.jw_v[end] + ssub.jw_ET[end]
    ssub.jw[end] = limit_lower_flux(water, jw*sbot.dt, θw, θwi, θsat, sat, Δz)/sbot.dt
end

function CryoGrid.boundaryflux(::Dirichlet, bc::WaterBC, top::Top, water::WaterBalance, sub::SubSurface, stop, ssub)
    Δk = CryoGrid.thickness(sub, ssub, first) # using `thickness` allows for generic layer implementations
    @inbounds let ψupper=boundaryvalue(bc, stop),
        ψsub=ssub.ψ[1],
        k=ssub.kw[1],
        δ=Δk/2; # distance to boundary
        Numerics.flux(ψupper, ψsub, δ, k)
    end
end
function CryoGrid.boundaryflux(::Dirichlet, bc::WaterBC, bot::Bottom, water::WaterBalance, sub::SubSurface, sbot, ssub)
    Δk = CryoGrid.thickness(sub, ssub, last) # using `thickness` allows for generic layer implementations
    @inbounds let ψlower=boundaryvalue(bc, sbot),
        ψsub=ssub.ψ[end],
        k=ssub.kw[end],
        δ=Δk/2; # distance to boundary
        Numerics.flux(ψsub, ψlower, δ, k)
    end
end

function CryoGrid.interact!(top::Top, bc::WaterBC, sub::SubSurface, water::WaterBalance, stop, ssub)
    ssub.jw_v[1] += CryoGrid.boundaryflux(bc, top, water, sub, stop, ssub)*ssub.dt
    balancefluxes!(top, bc, sub, water, stop, ssub)
    return nothing
end
function CryoGrid.interact!(sub::SubSurface, water::WaterBalance, bot::Bottom, bc::WaterBC, ssub, sbot)
    ssub.jw_v[end] += CryoGrid.boundaryflux(bc, bot, water, sub, ssub, sbot)*ssub.dt
    balancefluxes!(sub, water, bot, bc, ssub, sbot)
    return nothing
end

function CryoGrid.interact!(sub1::SubSurface, water1::WaterBalance{<:BucketScheme}, sub2::SubSurface, water2::WaterBalance{<:BucketScheme}, state1, state2)
    interact_ET!(sub1, water1, sub2, water2, state1, state2)
    θw₁ = state1.θw[end]
    θw₂ = state2.θw[1]
    θwi₁ = state1.θwi[end]
    θwi₂ = state2.θwi[1]
    θsat₁ = state1.θsat[end]
    θsat₂ = state2.θsat[1]
    sat₁ = state1.sat[end]
    sat₂ = state2.sat[1]
    Δz₁ = CryoGrid.thickness(sub1, state1, last)
    Δz₂ = CryoGrid.thickness(sub2, state2, first)
    # take minimum water content from upper layer where water would drain from
    θmin₁ = minwater(sub1, water1, state1, lastindex(state1.θw))
    kwc₁ = state1.kwc[end]
    kwc₂ = state2.kwc[1]
    kw = state1.kw[end] = state2.kw[1] = min(kwc₁, kwc₂)
    jw_v = advectiveflux(θw₁, θmin₁, kw)
    jw_ET = state2.jw_ET[1]
    jw = (jw_v + jw_ET)*state1.dt
    # setting both jw[end] on the upper layer and jw[1] on the lower layer is redundant since they refer to the same
    # element of the same underlying state array, but it's nice for clarity
    state1.jw[end] = state2.jw[1] = balanceflux(water1, water2, jw, θw₁, θw₂, θwi₁, θwi₂, θsat₁, θsat₂, sat₁, sat₂, Δz₁, Δz₂)/state1.dt
    return nothing
end
