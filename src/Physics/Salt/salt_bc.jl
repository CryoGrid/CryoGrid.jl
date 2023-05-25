Base.@kwdef struct SaltGradient{TbenthicSalt, TsurfaceState} <: BoundaryProcess{SaltMassBalance}
    benthicSalt::TbenthicSalt
    surfaceState::TsurfaceState
end

CryoGrid.BCKind(::Type{<:SaltGradient}) = CryoGrid.Dirichlet()

function CryoGrid.interact!(top::Top, bc::SaltGradient, soil::Soil, salt::SaltMassBalance, stop, ssoil)
    #upper boundary
    δ = CryoGrid.thickness(soil, ssoil, 1)
    surfaceState = bc.surfaceState(stop.t)
    benthicSalt = bc.benthicSalt(stop.t)
    flux = if (surfaceState == 0) # if is inundated
        -ssoil.dₛ[1] * (ssoil.c[1] - benthicSalt) / abs(δ / 2)
    else # under glacier or subaerial
        0.0
    end
    ssoil.∂c∂t[1] = flux / δ
    #lower boundary should get own function, but is zero anyhow so gets this by default
end
