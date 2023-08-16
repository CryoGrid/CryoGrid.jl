SaltProperties(
    τ = Param(1.5), # Turtuosity
    dₛ₀ = Param(8.0e-10, units=u"m^2/s"), # salt diffusion coefficient    
) = (; τ, dₛ₀)

abstract type SaltOperator end
# just one implementation for now
struct SaltDiffusion <: SaltOperator end

Base.@kwdef struct SaltMassBalance{TOp,Tdt,Tprop} <: SubSurfaceProcess
    op::TOp = SaltDiffusion()
    prop::Tprop = SaltProperties() # salt diffusion properties
    dtlim::Tdt = CryoGrid.CFL(0.5, CryoGrid.MaxDelta(5.0e-3)) # timestep limiter
end

Base.@kwdef struct SalineSoil{Tpara,Theat,Tsalt,Twater,Taux} <: AbstractGround{Tpara,Theat,Twater}
    para::Tpara = MineralOrganic()
    heat::Theat = HeatBalance(:T, freezecurve=DallAmicoSalt())
    salt::Tsalt = SaltMassBalance()
    water::Twater = WaterBalance(NoFlow())
    aux::Taux = nothing
end
SalineSoil(para::SoilParameterization; kwargs...) = SalineSoil(;para, kwargs...)

# type alias for coupled heat and salt diffusion
const CoupledHeatSalt{THeat,TSalt} = Coupled2{TSalt,THeat} where {THeat<:HeatBalance, TSalt<:SaltMassBalance}

# Soil methods
Soils.porosity(::SalineSoil, state) = state.por

# CryoGrid methods

CryoGrid.processes(soil::SalineSoil) = Coupled(soil.water, Coupled(soil.salt, soil.heat))

CryoGrid.variables(soil::SalineSoil) = (
    Diagnostic(:por, OnGrid(Cells), domain=0..1),
    CryoGrid.variables(soil, processes(soil))...,
)

CryoGrid.initializers(soil::SalineSoil) = (
    # default initializer for porosity and saturation
    initializer(:por, soil.para.por),
    initializer(:sat, soil.para.sat),
)
