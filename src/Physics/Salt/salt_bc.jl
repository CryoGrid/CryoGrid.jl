# type aliases for salt boundary conditions
const SaltBC = BoundaryProcess{T} where {SaltMassBalance<:T<:SubSurfaceProcess}
const SaltHeatBC{TSalt,THeat} = Coupled2{THeat,TSalt} where {TSalt<:SaltBC,THeat<:HeatBC}

SaltHeatBC(saltbc::SaltBC, heatbc::HeatBC) = Coupled(saltbc, heatbc)

Base.@kwdef struct SaltGradient{TbenthicSalt, TsurfaceState} <: BoundaryProcess{SaltMassBalance}
    benthicSalt::TbenthicSalt
    surfaceState::TsurfaceState
end

benthicSalt(bc::SaltGradient) = bc.benthicSalt

surfaceState(bc::SaltGradient) = bc.surfaceState

CryoGrid.BCKind(::Type{<:SaltGradient}) = CryoGrid.Dirichlet()

function CryoGrid.interact!(::Top, bc::SaltGradient, soil::SalineGround, ::SaltMassBalance, stop, ssed)
    # upper boundary
    surfaceState_t = surfaceState(bc)
    benthicSalt_t = benthicSalt(bc)
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
