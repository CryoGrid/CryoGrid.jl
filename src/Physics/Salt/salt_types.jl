SaltProperties(
    τ = Param(1.5), # Turtuosity
    dₛ₀ = Param(8.0e-10, units=u"m^2/s"), # salt diffusion coefficient    
) = (; τ, dₛ₀)


Base.@kwdef struct SaltMassBalance{TProp, Tdt} <: SubSurfaceProcess
    prop::TProp = SaltProperties() # salt diffusion properties
    dtlim::Tdt = CryoGrid.CFL(0.5, CryoGrid.MaxDelta(5.0e-3)) # timestep limiter
end

Base.@kwdef struct MarineSediment{Tpara,Theat,Tsalt,Twater,Tsp} <: Soil{Tpara,Theat,Twater}
    para::Tpara = MineralOrganic()
    heat::Theat = HeatBalance()
    salt::Tsalt = SaltMassBalance()
    water::Twater = nothing
    sp::Tsp = nothing
end
MarineSediment(para::SoilParameterization; kwargs...) = MarineSediment(para; kwargs...)

# Soil methods
Soils.porosity(::MarineSediment, state) = state.por

# CryoGrid methods
CryoGrid.variables(soil::MarineSediment) = (
    Diagnostic(:por, OnGrid(Cells), domain=0..1),
    CryoGrid.variables(soil, processes(soil))...,
)

CryoGrid.initializers(soil::MarineSediment) = (
    # default initializer that sets porosity state to the layer value
    initializer(:por, soil.para.por),
)
