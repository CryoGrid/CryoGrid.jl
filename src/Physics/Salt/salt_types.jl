SaltProperties(
    τ = param(1.5), # Turtuosity
    dₛ₀ = param(8.0e-10, units=u"m^2/s"), # salt diffusion coefficient    
) = (; τ, dₛ₀)

abstract type SaltOperator end
# just one implementation for now
struct SaltDiffusion <: SaltOperator end

Base.@kwdef struct SaltMassBalance{TOp,Tdt,Tprop} <: SubSurfaceProcess
    op::TOp = SaltDiffusion()
    prop::Tprop = SaltProperties() # salt diffusion properties
    dtlim::Tdt = CryoGrid.CFL(0.5, CryoGrid.MaxDelta(5.0e-3)) # timestep limiter
end

Base.@kwdef struct SalineGround{Tpara,Theat,Tsalt,Twater,Taux} <: AbstractGround{Tpara,Theat,Twater}
    para::Tpara = SimpleSoil(freezecurve=DallAmicoSalt())
    heat::Theat = HeatBalance(:T)
    salt::Tsalt = SaltMassBalance()
    water::Twater = WaterBalance(NoFlow())
    aux::Taux = nothing
end
SalineGround(para::SoilParameterization; kwargs...) = SalineGround(;para, kwargs...)

# type alias for coupled heat and salt diffusion
const CoupledHeatSalt{THeat,TSalt} = Coupled2{TSalt,THeat} where {THeat<:HeatBalance, TSalt<:SaltMassBalance}

# Soil methods
Soils.porosity(::SalineGround, state) = state.por

# CryoGrid methods

CryoGrid.processes(soil::SalineGround) = Coupled(soil.water, Coupled(soil.salt, soil.heat))

CryoGrid.variables(soil::SalineGround) = (
    Diagnostic(:por, OnGrid(Cells), domain=0..1),
    CryoGrid.variables(soil, processes(soil))...,
)

CryoGrid.initializers(soil::SalineGround) = (
    # default initializer for porosity and saturation
    initializer(:por, soil.para.por),
    initializer(:sat, soil.para.sat),
)
