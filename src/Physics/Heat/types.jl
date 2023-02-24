"""
    HeatOperator{progvar}

Base type for different numerical formulations of heat conduction.
"""
abstract type HeatOperator{progvar} end
"""
Type alias for `HeatOperator{:H}`, i.e. enthalpy-based heat conduction operators.
"""
const Enthalpy = HeatOperator{:H}
"""
Type alias for `HeatOperator{:T}`, i.e. temperature-based heat conduction operators.
"""
const Temperature = HeatOperator{:T}

thermalconductivity(op::HeatOperator) = op.cond
heatcapacity(op::HeatOperator) = op.hc

# Heat Balance type
"""
    HeatBalance{Tfc<:FreezeCurve,THeatOp<:HeatOperator,Tdt,Tprop} <: SubSurfaceProcess

Represents subsurface heat transfer processes. The formulation of heat transfer is governed by
the `HeatOperator`, `op` and 
"""
struct HeatBalance{Tfc<:FreezeCurve,THeatOp<:HeatOperator,Tdt,Tprop} <: SubSurfaceProcess
    op::THeatOp
    prop::Tprop
    freezecurve::Tfc
    dtlim::Tdt  # timestep limiter
end
# Heat operators
"""
    Diffusion{progvar,Tcond,Thc} <: HeatOperator

Represents a standard method-of-lines (MOL) forward diffusion operator for heat conduction with prognostic
`progvar`, typically either temperature `:T` or enthalpy (internal energy) `:H`.
"""
struct Diffusion{progvar,Tcond,Thc} <: HeatOperator{progvar}
    cond::Tcond
    hc::Thc
    Diffusion(progvar::Symbol, cond=quadratic_parallel_conductivity, hc=weighted_average_heatcapacity) = new{progvar,typeof(cond),typeof(hc)}(cond, hc)
end
# define constructor to allow for automatic reconstruction
ConstructionBase.constructorof(::Type{<:Diffusion{progvar}}) where {progvar} = (cond,hc) -> Diffusion(progvar, cond, hc)
"""
    EnthalpyImplicit <: HeatOperator{:H}

Implicit enthalpy formulation of Swaminathan and Voller (1992) and Langer et al. (2022). Note that this
heat operator formulation does not compute a divergence `∂H∂t` but only computes the necessary diffusion
coefficients for use by an appropriate solver. See the `Solvers.LiteImplicit` module for the appropriate
solver algorithms.
"""
struct EnthalpyImplicit{Tcond,Thc} <: HeatOperator{:H}
    cond::Tcond
    hc::Thc
    EnthalpyImplicit(cond=quadratic_parallel_conductivity, hc=weighted_average_heatcapacity) = new{typeof(cond),typeof(hc)}(cond, hc)
end

"""
Numerical constants for pararameterizing heat processes.
"""
Utils.@properties HeatBalanceProperties(
    ρw = CryoGrid.Constants.ρw,
    Lsl = CryoGrid.Constants.Lsl,
    L = ρw*Lsl,
)
# do not parameterize heat properties
CryoGrid.parameterize(prop::HeatBalanceProperties) = prop

"""
    ThermalProperties

Basic material thermal properties.
"""
Utils.@properties ThermalProperties(
    kh_w = 0.57u"J/s/m/K", # thermal conductivity of water [Hillel (1982)]
    kh_i = 2.2u"J/s/m/K", # thermal conductivity of ice [Hillel (1982)]
    kh_a = 0.025u"J/s/m/K", # thermal conductivity of air [Hillel (1982)]
    ch_w = 4.2e6u"J/K/m^3", # heat capacity of water
    ch_i = 1.9e6u"J/K/m^3", # heat capacity of ice
    ch_a = 0.00125e6u"J/K/m^3", # heat capacity of air
)
function CryoGrid.parameterize(props::ThermalProperties)
    @set! props.kh_w = CryoGrid.parameterize(props.kh_w)
    @set! props.kh_i = CryoGrid.parameterize(props.kh_i)
    @set! props.kh_a = CryoGrid.parameterize(props.kh_a)
    @set! props.ch_w = CryoGrid.parameterize(props.ch_w)
    @set! props.ch_i = CryoGrid.parameterize(props.ch_i)
    @set! props.ch_a = CryoGrid.parameterize(props.ch_a)
    return props
end

# default step limiters
default_dtlim(::Temperature) = CryoGrid.CFL(maxdelta=CryoGrid.MaxDelta(Inf))
default_dtlim(::Enthalpy) = CryoGrid.MaxDelta(1u"MJ")
default_dtlim(::HeatOperator) = nothing

# convenience constructors for specifying prognostic variable as symbol
HeatBalance(var::Symbol=:H; kwargs...) = HeatBalance(Val{var}(); kwargs...)
HeatBalance(::Val{:H}; kwargs...) = HeatBalance(Diffusion(:H); kwargs...)
HeatBalance(::Val{:T}; kwargs...) = HeatBalance(Diffusion(:T); kwargs...)
HeatBalance(op; freezecurve=FreeWater(), prop=HeatBalanceProperties(), dtlim=default_dtlim(op)) = HeatBalance(op, prop, deepcopy(freezecurve), dtlim)
HeatBalance(op::Temperature; freezecurve, prop=HeatBalanceProperties(), dtlim=default_dtlim(op)) = HeatBalance(op, prop, deepcopy(freezecurve), dtlim)
