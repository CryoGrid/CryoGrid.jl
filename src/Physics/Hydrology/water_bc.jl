# Boundary condition type aliases
const Dirichlet = CryoGrid.Dirichlet
const Neumann = CryoGrid.Neumann
const WaterBC = BoundaryProcess{T} where {WaterBalance<:T<:SubSurfaceProcess}

# Constant boundary constructors
ConstantInfiltration(value::Quantity) = ConstantBC(WaterBalance, Neumann, uconvert(u"m/s", value))
ConstantInfiltration(value) = ConstantBC(WaterBalance, Neumann, value)
ImpermeableBoundary() = ConstantBC(WaterBalance, Neumann, 0.0u"m/s")

"""
    Rainfall{Train<:Forcing{u"m/s"}} <: BoundaryProcess{WaterBalance}

Basic rainfall boundary condition for `WaterBalance` which simply invokes the given precipitation
forcing at the current time `t`.
"""
struct Rainfall{Train<:Forcing{u"m/s"}} <: BoundaryProcess{WaterBalance}
    rain::Train
end

CryoGrid.BCKind(::Type{<:Rainfall}) = Neumann()

function CryoGrid.boundaryvalue(bc::Rainfall, ::Top, ::WaterBalance, sub::SubSurface, stop, ssub)
    rainfall_rate = bc.rain(stop.t)
    # take the minimum of the current rainfall rate and the hydraulic conductivity at the top of the upper grid cell;
    # note that this assumes rainfall to be in m/s
    return min(rainfall_rate, ssub.kw[1])
end

@inline CryoGrid.boundaryflux(::Neumann, bc::WaterBC, top::Top, water::WaterBalance, sub::SubSurface, stop, ssub) = boundaryvalue(bc,top,water,sub,stop,ssub)
@inline CryoGrid.boundaryflux(::Neumann, bc::WaterBC, bot::Bottom, water::WaterBalance, sub::SubSurface, sbot, ssub) = boundaryvalue(bc,bot,water,sub,sbot,ssub)
@inline function CryoGrid.boundaryflux(::Dirichlet, bc::WaterBC, top::Top, water::WaterBalance, sub::SubSurface, stop, ssub)
    Δk = CryoGrid.thickness(sub, ssub, first) # using `thickness` allows for generic layer implementations
    @inbounds let ψupper=boundaryvalue(bc,top,water,sub,stop,ssub),
        ψsub=ssub.ψ[1],
        k=ssub.kw[1],
        δ=Δk/2; # distance to boundary
        -k*(ψsub-ψupper)/δ
    end
end
@inline function CryoGrid.boundaryflux(::Dirichlet, bc::WaterBC, bot::Bottom, water::WaterBalance, sub::SubSurface, sbot, ssub)
    Δk = CryoGrid.thickness(sub, ssub, last) # using `thickness` allows for generic layer implementations
    @inbounds let ψlower=boundaryvalue(bc,bot,water,sub,sbot,ssub),
        ψsub=ssub.ψ[end],
        k=ssub.kw[end],
        δ=Δk/2; # distance to boundary
        # note again the inverted sign; positive here means *upward from* the bottom boundary
        # TODO: maybe change this convention? it seems needlessly confusing and bug-prone.
        k*(ψlower-ψsub)/δ
    end
end

function CryoGrid.interact!(top::Top, bc::WaterBC, sub::SubSurface, water::WaterBalance, stop, ssub)
    θw = ssub.θw[1]
    θwi = ssub.θwi[1]
    θsat = ssub.θsat[1]
    sat = ssub.sat[1]
    Δz = CryoGrid.thickness(sub, ssub, first)
    jw_up = CryoGrid.boundaryflux(bc, top, water, sub, stop, ssub)*ssub.dt
    jw_up = limit_upper_flux(water, jw_up, θw, θwi, θsat, sat, Δz)
    ssub.jw[1] += jw_up
    return nothing
end
function CryoGrid.interact!(sub::SubSurface, water::WaterBalance, bot::Bottom, bc::WaterBC, ssub, sbot)
    θw = ssub.θw[end]
    θwi = ssub.θwi[end]
    θsat = ssub.θsat[end]
    sat = ssub.sat[end]
    Δz = CryoGrid.thickness(sub, ssub, first)
    # sign flipped due to positive downward convention
    jw_lo = -CryoGrid.boundaryflux(bc, bot, water, sub, ssub, sbot)*ssub.dt
    jw_lo = limit_lower_flux(water, jw_lo, θw, θwi, θsat, sat, Δz)
    ssub.jw[end] += jw_lo
    return nothing
end