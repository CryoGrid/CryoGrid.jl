# type alias for salt boundary conditions
const SaltBC = BoundaryProcess{T} where {SaltMassBalance<:T<:SubSurfaceProcess}

Base.@kwdef struct SaltGradient{TbenthicSalt, TsurfaceState} <: BoundaryProcess{SaltMassBalance}
    benthicSalt::TbenthicSalt
    surfaceState::TsurfaceState
end

benthicSalt(bc::SaltGradient, t) = bc.benthicSalt
benthicSalt(bc::SaltGradient{<:Forcing}, t) = bc.benthicSalt(t)

surfaceState(bc::SaltGradient, t) = bc.surfaceState
surfaceState(bc::SaltGradient{<:Forcing}, t) = bc.surfaceState(t)

CryoGrid.BCKind(::Type{<:SaltGradient}) = CryoGrid.Dirichlet()

function CryoGrid.interact!(::Top, bc::SaltGradient, soil::SaltySoil, ::SaltMassBalance, stop, ssed)
    # upper boundary
    surfaceState_t = surfaceState(bc, stop.t)
    benthicSalt_t = benthicSalt(bc, stop.t)
    # distance to boundary; 1/2 thickness of first grid cell
    dz = CryoGrid.thickness(soil, ssed, first) / 2
    flux = if (surfaceState_t == 0) # if is inundated
        -ssed.dₛ[1] * (ssed.c[1] - benthicSalt_t) / dz
    else # under glacier or subaerial
        0.0
    end
    ssed.jc[1] = flux
    #lower boundary should get own function, but is zero anyhow so gets this by default
end

Base.@kwdef struct HeatSaltBC{Theat,Tsalt} <: BoundaryProcess{Union{SaltMassBalance,HeatBalance}}
    heat::Theat
    salt::Tsalt
end

function CryoGrid.updatestate!(top::Top, bc::HeatSaltBC, state)
    updatestate!(top, bc.salt, state)
    updatestate!(top, bc.heat, state)
end

function CryoGrid.interact!(top::Top, bc::HeatSaltBC, sub::SubSurface, ps::Coupled(SaltMassBalance, HeatBalance), stop, ssub)
    interact!(top, bc.salt, sub, ps[1], stop, ssub)
    interact!(top, bc.heat, sub, ps[2], stop, ssub)
end

function CryoGrid.computefluxes!(top::Top, bc::HeatSaltBC, state)
    computefluxes!(top, bc.salt, state)
    computefluxes!(top, bc.heat, state)
end

