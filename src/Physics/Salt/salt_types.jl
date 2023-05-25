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

Base.@kwdef struct MarineSediment{Tpara,Theat,Tsalt,Twater,Tsp} <: Soil{Tpara,Theat,Twater}
    para::Tpara = MineralOrganic()
    heat::Theat = HeatBalance(:T, freezecurve=DallAmicoSalt())
    salt::Tsalt = SaltMassBalance()
    water::Twater = nothing
    sp::Tsp = nothing
end
MarineSediment(para::SoilParameterization; kwargs...) = MarineSediment(para; kwargs...)

# Soil methods
Soils.porosity(::MarineSediment, state) = state.por

# CryoGrid methods

CryoGrid.processes(sediment::MarineSediment) = Coupled(sediment.salt, sediment.heat)
CryoGrid.processes(sediment::MarineSediment{Tpara,Theat,Tsalt,<:WaterBalance}) where {Tpara,Theat,Tsalt} = Coupled(sediment.water, sediment.salt, sediment.heat)

CryoGrid.variables(sediment::MarineSediment) = (
    Diagnostic(:por, OnGrid(Cells), domain=0..1),
    CryoGrid.variables(sediment, processes(sediment))...,
)

CryoGrid.initializers(sediment::MarineSediment) = (
    # default initializer that sets porosity state to the layer value
    initializer(:por, sediment.para.por),
)
