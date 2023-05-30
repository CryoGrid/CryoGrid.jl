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

function CryoGrid.interact!(::Top, bc::SaltGradient, soil::Soil, ::SaltMassBalance, stop, ssoil)
    # upper boundary
    surfaceState = surfaceState(bc, stop.t)
    benthicSalt = benthicSalt(bc, stop.t)
    flux = if (surfaceState == 0) # if is inundated
        -ssoil.dₛ[1] * (ssoil.c[1] - benthicSalt) / abs(δ / 2)
    else # under glacier or subaerial
        0.0
    end
    ssoil.jc[1] = flux
    #lower boundary should get own function, but is zero anyhow so gets this by default
end
